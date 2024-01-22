#' Function to extract positions of methylated cytosines from the BAM file in forward reference coordinates.
#'
#' @param bamfile A path to a BAM file containing methylation specific tags (MM and ML).
#' @param bamindex A path to the BAM file index. By default index file is expected to be present in the same directory. If not present it will be created.
#' @param region A user defined region from where the methylation data will be extracted.
#' @param min.mapq A minimum mapping quality of aligned reads to be extracted from provided BAM file.
#' @param add.phase Set to \code{TRUE} if phasing specific tag 'HP' should be added if defined in the BAM file.
#' @importFrom Rsamtools indexBam ScanBamParam BamFile scanBam
#' @importFrom methods is as
#' @importFrom GenomicAlignments readGAlignments mapFromAlignments
#' @importFrom GenomicRanges sort disjointBins GRanges strand start end
#' @importFrom GenomeInfoDb seqnames
#' @importFrom IRanges IRanges reflect subsetByOverlaps ranges
#' @importFrom Biostrings reverseComplement matchPattern
#' @return A \code{list} of \code{\link{GRanges-class}} objects called `reads` storing read boundaries and positions
#' and `methylCs` reporting positions of methylated cytosines in reference coordinates.
#' @author David Porubsky
#' @export
#' @examples
#' ## Get BAM to extract methylation from
#' bam.file <- system.file("extdata", "test.bam", package = "ggNucFreak")
#' ## Set region to extract methylation from
#' region <- as('chr11:2787889-2790647', "GRanges")
#' region <- GenomicRanges::resize(region, width = GenomicRanges::width(region) + (GenomicRanges::width(region) * 2), fix = 'center')
#' ## Extract methylation
#' mCs <- extractMethylCs(bamfile = bam.file, region=region)
#' mCs
#'
extractMethylCs <- function(bamfile=NULL, bamindex=bamfile, region=NULL, min.mapq=10, add.phase=FALSE) {
  ## Check user input
  bamindex.raw <- sub('\\.bai$', '', bamindex)
  bamindex <- paste0(bamindex.raw,'.bai')
  if (!file.exists(bamindex)) {
    bamindex.own <- Rsamtools::indexBam(bamfile)
    warning("Couldn't find BAM index-file ", bamindex, ". Indexing ...")
    bamindex <- bamindex.own
  }

  if (!is.null(region)) {
    if (is.character(region)) {
      region.gr <- methods::as(region, "GRanges")
    } else if (is(region, "GRanges")) {
      region.gr <- region
    } else {
      region.gr <- NULL
      message("Parameter 'region' can either be 'GRanges' object or character string 'chr#:start-end'!!!")
    }
  }

  ## Read in BAM [Genomic Alignments]
  if (!is.null(region)) {
    aln.data <- GenomicAlignments::readGAlignments(bamfile, index=bamindex,
                                                   param=Rsamtools::ScanBamParam(which=range(region.gr),
                                                                                 what=c('mapq', 'flag', 'seq', 'qname')))
  } else {
    aln.data <- GenomicAlignments::readGAlignments(bamfile, index=bamindex,
                                                   param=Rsamtools::ScanBamParam(what=c('mapq', 'flag', 'seq', 'qname')))
  }
  ## Sort alignments
  aln.data <- GenomicRanges::sort(aln.data, ignore.strand = TRUE)

  ## Define BAM tags to load
  if (add.phase) {
    tags <- c('MM', 'ML', 'HP')
  } else {
    tags <- c('MM', 'ML')
  }
  ## Load methylation specific tags and other user defined tags [Rsamtools]
  if (!is.null(region)) {
    params <- Rsamtools::ScanBamParam(tag=tags, which=range(region.gr))
  } else {
    params <- Rsamtools::ScanBamParam(tag=tags)
  }
  bamFile <- Rsamtools::BamFile(bamfile)
  bam.tags <- Rsamtools::scanBam(bamfile, param = params)

  ## Add methylation tags
  data.gr <- methods::as(aln.data, 'GRanges')
  data.gr$MM <- bam.tags[[1]]$tag$MM
  data.gr$ML <- sapply(bam.tags[[1]]$tag$ML, function(x) paste(x, collapse = ','))

  ## Add phasing tag
  if (add.phase) {
    if (methods::is(bam.tags[[1]]$tag$HP, 'numeric')) {
      data.gr$HP <-  bam.tags[[1]]$tag$HP
    } else {
      add.phase <- FALSE
      warning("Haplotype specific BAM tag 'HP' is likely not defined in submitted bam !!!")
    }
  }

  ## Keep only reads with reported methylation
  mask.methyl <- !is.na(data.gr$MM) & data.gr$MM != ''
  data.gr <- data.gr[mask.methyl]
  aln.data <- aln.data[mask.methyl]

  ## Get plotting coordinates
  ## Get level of all overlapping reads
  data.gr$level <- GenomicRanges::disjointBins(data.gr, ignore.strand = TRUE)
  ## Get level of overlapping reads by phase
  if (add.phase) {
    data.gr$phase.level <- NA
    data.gr[data.gr$HP == 1]$phase.level <- GenomicRanges::disjointBins(data.gr[data.gr$HP == 1], ignore.strand = TRUE)
    data.gr[data.gr$HP == 2]$phase.level <- GenomicRanges::disjointBins(data.gr[data.gr$HP == 2], ignore.strand = TRUE)
  }

  all.reads <- list()
  for (i in seq_along(data.gr)) {
    ## Process single read
    read.gr <- data.gr[i]
    aln <- aln.data[i]
    ## Get read sequence
    seq <- unlist(read.gr$seq)
    if (as.character(GenomicRanges::strand(read.gr)) == '-') {
      seq <- Biostrings::reverseComplement(seq)
    }

    ## Detect positions of C nucleotides
    c.pos <- Biostrings::matchPattern('C', seq)
    ## Detect positions of CG dinucleotides
    cg.pos <- Biostrings::matchPattern('CG', seq)

    ## Get position of methylated cytosines
    ## Get MM tag ##
    mm.tag <- strsplit(read.gr$MM, ',|;')[[1]]
    mm.tag <- as.numeric(mm.tag[-1])
    # if (as.character(strand(read.gr)) == '-') { ## Not needed as modified base is shown in original read orientation
    #   mm.tag <- rev(mm.tag)
    # }
    modifs <- rep(1, length(mm.tag))
    ## Get cumulative sum of unmethylated C's and following methylated C's
    cumsum.pos <- cumsum( c(rbind(mm.tag, modifs)) )
    ## Get methylated C position by taking even positions in cumsum.pos
    methC.pos <- cumsum.pos[1:length(cumsum.pos) %% 2 == 0]
    ## Create GRanges object of cytosine modification
    c.pos.gr <- GenomicRanges::GRanges(seqnames = GenomeInfoDb::seqnames(read.gr), ranges = c.pos@ranges, qname = read.gr$qname, strand = GenomicRanges::strand(read.gr))
    c.pos.gr$label <- 'C'
    c.pos.gr$label[methC.pos] <- 'mC'
    ## Mark CpG positions
    c.pos.gr$CpG <- FALSE
    c.pos.gr$CpG[start(c.pos.gr) %in% start(cg.pos)] <- TRUE

    ## Get ML tag ##
    ml.tag <- strsplit(read.gr$ML, ',|;')[[1]]
    ml.tag <- as.numeric(ml.tag)
    ## Convert to probabilities
    lower.prob <- (ml.tag / 256)
    upper.prob <- ((ml.tag + 1) / 256)
    ml.prob.tag <- lower.prob + ((upper.prob - lower.prob) / 2)
    c.pos.gr$m.prob <- NA
    c.pos.gr$m.prob[methC.pos] <- ml.prob.tag

    ## Add phasing tag
    if (add.phase) {
      c.pos.gr$HP <- read.gr$HP
    }

    ## Flip the orientation to the reference (forward)
    if (as.character(GenomicRanges::strand(read.gr)) == '-') {
      bounds <- IRanges::IRanges(start = 0L, end = length(seq))
      IRanges::ranges(c.pos.gr) <- IRanges::reflect(IRanges::ranges(c.pos.gr), bounds = bounds)
    }

    ## Map cytosine positions to alignment
    alignments <- rep(aln, length(c.pos.gr))
    names(alignments) <- 1:length(c.pos.gr)
    names(c.pos.gr) <- 1:length(c.pos.gr)
    c.aln.gr <- GenomicAlignments::mapFromAlignments(x = c.pos.gr, alignments = alignments)
    IRanges::ranges(c.pos.gr) <- IRanges::ranges(c.aln.gr)
    names(c.pos.gr) <- NULL
    #
    # bounds <- IRanges::IRanges(start = start(aln), end = end(aln))
    # ranges(c.pos.gr) <- reflect(ranges(c.pos.gr), bounds = bounds)

    ## Sort by position
    c.pos.gr <- GenomicRanges::sort(c.pos.gr)

    ## For testing
    #subsetByOverlaps(c.pos.gr, region)

    ## Add plotting level
    c.pos.gr$level <- read.gr$level
    if (add.phase) {
      c.pos.gr$phase.level <- read.gr$phase.level
    }

    ## Store ranges
    all.reads[[i]] <- c.pos.gr
  }
  methylCs.gr <- do.call(c,  all.reads)

  ## Subset data to a user specified region
  if (methods::is(region.gr, 'GRanges')) {
    methylCs.gr <- IRanges::subsetByOverlaps(methylCs.gr, region.gr)
    data.gr$start.crop <- pmax(GenomicRanges::start(data.gr), GenomicRanges::start(region.gr))
    data.gr$end.crop <- pmin(GenomicRanges::end(data.gr), GenomicRanges::end(region.gr))
  }

  ## Return final object
  names(methylCs.gr) <- NULL
  names(data.gr) <- NULL
  return(list(reads=data.gr, methylCs=methylCs.gr))
}
