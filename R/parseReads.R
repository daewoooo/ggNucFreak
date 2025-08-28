#' Function to extract and parse reads from the BAM file with respect to the reference.
#'
#' @param min.insertion.size A minimum size (in base pairs) of an insertion to be retained.
#' @param min.deletion.size A minimum size (in base pairs) of a deletion to be retained.
#' @param collapse.mismatches Set to \code{TRUE} if mismatches should be collapsed in order expand matched regions.
#' @inheritParams exportMismatches
#' @return A \code{list} containing two GRanges objects called `reads` storing reported mismatches
#' and `alignments` reporting overall coverage per genomic position.
#' @author David Porubsky
#' @export
#' 
parseReads <- function(bamfile=NULL, bamindex=bamfile, secondary.aln=FALSE, keep.duplicates=FALSE, region=NULL, min.mapq=10, min.baseq=20, by.phase=FALSE, min.deletion.size = 10, min.insertion.size = 10, collapse.mismatches = TRUE) {
  ## Check user input ##
  ## Load BAM file index
  bamindex.raw <- sub('\\.bai$', '', bamindex)
  bamindex <- paste0(bamindex.raw,'.bai')
  if (!file.exists(bamindex)) {
    bamindex.own <- Rsamtools::indexBam(bamfile)
    warning("Couldn't find BAM index-file ", bamindex, ". Indexing ...")
    bamindex <- bamindex.own
  }
  
  ## Load user defined region of interest
  if (!is.null(region)) {
    if (is.character(region)) {
      region.gr <- methods::as(region, "GRanges")
    } else if (methods::is(region, "GRanges")) {
      region.gr <- region
    } else {
      stop("Parameter 'region' MUST be defined as either 'GRanges' object or character string 'chr#:start-end'!!!")
    }
  } else {
    stop("Parameter 'region' MUST be defined as either 'GRanges' object or character string 'chr#:start-end'!!!")
  }
  ## Restrict seqlevels
  region.gr <- GenomeInfoDb::keepSeqlevels(region.gr, value = as.character(GenomeInfoDb::seqnames(region.gr)), pruning.mode = 'coarse')
  
  ## Read in BAM [Genomic Alignments] ##
  ######################################
  ## Define BAM flags to filter
  flags <- Rsamtools::scanBamFlag(isSecondaryAlignment = secondary.aln,
                                  isDuplicate = keep.duplicates)
  
  if (by.phase) {
    aln.data <- GenomicAlignments::readGAlignments(bamfile, index=bamindex,
                                                   param=Rsamtools::ScanBamParam(which = region.gr,
                                                                                 flag = flags,
                                                                                 what=c('mapq', 'flag', 'seq', 'qual', 'qname'),
                                                                                 tag = 'HP'))
    ## Keep only alignments with defined haplotype tag
    aln.data <- aln.data[!is.na(mcols(aln.data)$HP)]
  } else {
    aln.data <- GenomicAlignments::readGAlignments(bamfile, index=bamindex,
                                                   param=Rsamtools::ScanBamParam(which = region.gr,
                                                                                 flag = flags,
                                                                                 what=c('mapq', 'flag', 'seq', 'qual', 'qname')))
    mcols(aln.data)$HP <- 0
  }
  
  ## Convert to GRanges object
  mcols(aln.data)$cigar <- aln.data@cigar
  aln.data <- as(aln.data, 'GRanges')
  
  ## Filter by mapping quality
  if (min.mapq > 0) {
    aln.data <- aln.data[aln.data$mapq >= min.mapq]
  }
  ## Sort alignments
  aln.data <- GenomicRanges::sort(aln.data, ignore.strand = TRUE)
  ## Restrict seqlevels to region of interest
  aln.data <- GenomeInfoDb::keepSeqlevels(aln.data, value = as.character(GenomeInfoDb::seqnames(region.gr)), pruning.mode = 'coarse')
  
  ## Parse alignments per haplotype
  haps <- unique(aln.data$HP)
  aln.ranges.l <- list()
  reads.l <- list()
  for (i in seq_along(haps)) {
    hap <- haps[i]
    sub.data <- aln.data[aln.data$HP == hap]
    sub.data <- sort(sub.data, ignore.strand=TRUE)
    seqnames <- GenomeInfoDb::seqlevels(sub.data)
    ## Get plotting level for each read
    sub.data$level <- disjointBins(sub.data, ignore.strand=TRUE)
    
    ## Parse CIGAR string of each read
    aln.ranges <- list()
    for (j in seq_along(sub.data)) {
      aln <- sub.data[j,]
      qname <- aln$qname
      
      ## Get matched bases on reference space
      aln.ref <- SVbyEye::parseCigarString(cigar.str = aln$cigar, coordinate.space = 'reference')
      match <- aln.ref$match
      mismatch <- aln.ref$mismatch
      deletion <- aln.ref$deletion
      ## Matched ranges
      if (!is.null(match)) {
        match.gr <- GRanges(seqnames = seqnames(aln), ranges = match)
      } else {
        match.gr <- GRanges()
      }  
      ## Mismatched ranges
      if (!is.null(mismatch)) {
        ## Make sure all are one base changes
        pos <- sort(unique(c(start(mismatch), end(mismatch))))
        mismatch.gr <- GRanges(seqnames = seqnames(aln), ranges = IRanges(start = pos, end = pos))
      } else {
        mismatch.gr <- GRanges()
      } 
      ## Deleted ranges (missing in the query/reads)
      if (!is.null(deletion)) {
        del.gr <- GRanges(seqnames = seqnames(aln), ranges = deletion)
        del.mask <- GenomicRanges::width(del.gr) >= min.deletion.size
        del2reduce <- del.gr[!del.mask]
        del.gr <- del.gr[del.mask]
        if (length(del.gr) > 0) {
          del.gr$aln.type <- 'DEL'
          del.gr$aln.len <- width(del.gr)
          del.gr$seq <- NA
        } else {
          del.gr <- GRanges()
        }
      } else {
        del.gr <- GRanges()
        del2reduce <- GRanges()
      } 
      ## Collapse mismatched position if desired
      if (collapse.mismatches) {
        match.gr <- reduce(c(match.gr, mismatch.gr))
      } else {
        if (length(mismatch.gr) > 0) {
          piles <- GenomicAlignments::pileLettersAt(aln$seq, GenomeInfoDb::seqnames(aln), 1L, aln$cigar, mismatch.gr)
          mismatch.gr$aln.type <- 'mismatch'
          mismatch.gr$aln.len <- width(mismatch.gr)
          mismatch.gr$seq <- as.character(piles)
          ## Filter mismatches by base quality if defined
          if (min.baseq > 0) {
            quals <- GenomicAlignments::pileLettersAt(aln$qual, GenomeInfoDb::seqnames(aln), 1L, aln$cigar, mismatch.gr)
            quals <- sapply(as.character(quals), function(x) as.numeric(charToRaw(x))-33)
            mask.qual <- which(quals >= min.baseq)
            mismatch.gr <-  mismatch.gr[mask.qual]
          }
        }  
      }
      ## Collapse filtered deletions
      match.gr <- reduce(c(match.gr, del2reduce))
      match.gr$aln.type <- 'match'
      match.gr$aln.len <- width(match.gr)
      match.gr$seq <- NA
      
      ## Get inserted bases on query space
      aln.qry <- SVbyEye::parseCigarString(cigar.str = aln$cigar, coordinate.space = 'query')
      qry.insertion <- aln.qry$insertion
      ref.insertion <- aln.ref$insertion
      if (!is.null(qry.insertion)) {
        mask <- which(width(qry.insertion) > min.insertion.size)
        #ins.gr <- GRanges(seqnames = seqnames(aln), ranges = ref.insertion[mask])
        #ins.gr <- ins.gr[GenomicRanges::width(ins.gr) >= min.insertion.size]
        if (length(mask) > 0) {
          ins.gr <- GRanges(seqnames = seqnames(aln), ranges = ref.insertion[mask])
          ins.gr$aln.type <- 'INS'
          ins.gr$aln.len <- width(qry.insertion[mask])
          #ins.gr$aln.len <- NA
        } else {
          ins.gr <- GRanges()
        }
      } else {
        ins.gr <- GRanges()
      }
      
      ## Merge and sort ranges
      if (collapse.mismatches) {
        ranges <- c(match.gr, del.gr, ins.gr)
      } else {
        ranges <- c(match.gr, mismatch.gr, del.gr, ins.gr)
      }  
      ranges <- sort(ranges, ignore.strand=TRUE)
      ## Convert to reference coordinates
      ranges <- GenomicRanges::shift(ranges, shift = start(aln) - 1)
      ## Add extra info
      ranges$qname <- qname
      ranges$HP <- hap
      ## Add plotting level
      ranges$level <- aln$level
      ## Store
      aln.ranges[[j]] <- ranges
    }  
    aln.ranges.gr <- do.call(c, aln.ranges)
    
    ## Restrict ranges to region of interest
    aln.ranges.gr <- subsetByOverlaps(aln.ranges.gr, region.gr)
    start(aln.ranges.gr) <- pmax(start(aln.ranges.gr), start(region.gr))
    end(aln.ranges.gr) <- pmin(end(aln.ranges.gr), end(region.gr))
    start(sub.data) <- pmax(start(sub.data), start(region.gr))
    end(sub.data) <- pmin(end(sub.data), end(region.gr))
    
    ## For testing
    aln.df <- as.data.frame(aln.ranges.gr)
    reads.df <- as.data.frame(sub.data)
    match.df <- aln.df[aln.df$aln.type == 'match',]
    mismatch.df <- aln.df[aln.df$aln.type == 'mismatch',]
    ins.df <- aln.df[aln.df$aln.type == 'INS',]
    dels.df <- aln.df[aln.df$aln.type == 'DEL',]
    base.cols <- c('A' = 'coral2', 'C' = 'lightgreen', 'G' = 'steelblue3', 'T' = 'tan2')
    ggplot() +
      geom_segment(data = reads.df, aes(x=start, xend=end, y=level + 0.5, yend=level + 0.5)) +
      geom_rect(data = match.df, aes(xmin=start, xmax=end, ymin=level, ymax=level+1)) +
      geom_rect(data = reads.df, aes(xmin=start, xmax=end, ymin=level, ymax=level+1), color='white', fill=NA) +
      geom_rect(data = mismatch.df, aes(xmin=start, xmax=end, ymin=level, ymax=level+1, color=seq), fill=NA) +
      scale_color_manual(values = base.cols) +
      scale_y_continuous(expand = c(0, 0)) +
      facet_grid(HP ~ ., space='free', scale='free') +
      theme_minimal()
    
    ## Store read ranges per haplotype
    aln.ranges.l[[i]] <- aln.ranges.gr
    reads.l[[i]] <- sub.data[,names(mcols(sub.data)) %in% c('mapq', 'flag', 'qname', 'HP', 'level')]
  }
  all.reads <- do.call(c, reads.l)
  all.aln.ranges <- do.call(c, aln.ranges.l)
  ## Return ranges
  gr.list <- list(reads = all.reads, alignments = all.aln.ranges) 
  return(gr.list)
}  

## Function to be finished
plotReads <- function(reads.obj) {
  ## Check user input ##
  if (!methods::is(reads.obj, 'list')) {
    stop("Submitted object has to be a list containing genomic ranges 'reads' and 'alignments' exported by 'parseReads' function !!!")
    if (!all(names(reads.obj) %in% c("reads", "coverage"))) {
      stop("Submitted object does not contain expected genomic ranges 'reads' and 'alignments' exported by 'parseReads' function !!!")
    }
  }
  
  ## Get data to plot
  reads.df <- as.data.frame(reads.obj$reads)
  aln.df <- as.data.frame(reads.obj$alignments)
  
  ## Add haplotype labels
  if (!is.null(reads.df)) {
    reads.df$HP <- ifelse(reads.df$HP == 0, 'Unphased', paste0('Haplotype ', reads.df$HP))
  }
  if (!is.null(aln.df)) {
    aln.df$HP <- ifelse(aln.df$HP == 0, 'Unphased', paste0('Haplotype ', aln.df$HP))
  }
  
  ## Make a plot ##
  #################
  match.df <- aln.df[aln.df$aln.type == 'match',]
  ins.df <- aln.df[aln.df$aln.type == 'INS',]
  dels.df <- aln.df[aln.df$aln.type == 'DEL',]
  
  ## Set plottting theme
  custom_theme <- theme(panel.background = element_blank())
  
  plt <- ggplot() +
    geom_segment(data = reads.df, aes(x=start, xend=end, y=level + 0.5, yend=level + 0.5)) +
    geom_rect(data = match.df, aes(xmin=start, xmax=end, ymin=level, ymax=level+1)) +
    geom_rect(data = reads.df, aes(xmin=start, xmax=end, ymin=level, ymax=level+1), color='white', fill=NA) +
    scale_y_continuous(expand = c(0, 0), name = 'Read coverage') +
    scale_x_continuous(expand = c(0, 0), name = 'Genomic position (bp)') +
    facet_grid(HP ~ ., space='free', scale='free') +
    custom_theme
  
  # plt <- ggplot() +
  #   geom_segment(data = reads.df, aes(x=start, xend=end, y=level, yend=level)) +
  #   geom_segment(data = match.df, aes(x=start, xend=end, y=level, yend=level), size=3, lineend = 'round') +
  #   scale_y_continuous(expand = c(0.1, 0.1), name = 'Read coverage') +
  #   scale_x_continuous(expand = c(0, 0), name = 'Genomic position (bp)') +
  #   facet_grid(HP ~ ., space='free', scale='free') +
  #   custom_theme
  
  ## Add insertion positions
  if (nrow(ins.df) > 0) {
    plt <- plt +
      #geom_rect(data = reads.df, aes(xmin=start, xmax=end, ymin=level, ymax=level+1), color='white', fill=NA) +
      #geom_richtext(data = ins.df, aes(x=end, y=level + 0.5, label=SVLEN), inherit.aes = FALSE, size=1) +
      geom_label(data = ins.df, aes(x=end, y=level + 0.5, label=SVLEN), size=2)
  }
  
  ## Return final plot
  return(plt)
}
