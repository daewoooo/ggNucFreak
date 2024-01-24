#' Function to extract positions of mismatched bases from the BAM file with respect to the reference.
#'
#' @param bamfile A path to a BAM file.
#' @param add.phase Set to \code{TRUE} if phasing specific tag 'HP' should be added if defined in the BAM file.
#' @param min.baseq A minimum base quality per aligned base to be retained.
#' @param min.mismatches A minimum number of mismatch coverage to be reported.
#' @param by.phase If set to \code{TRUE} mismatches will be reported per haplotype specific BAM tag ('HP').
#' @param ref.fasta A path to a reference FASTA file to extract reference bases from.
#' @param bsgenome A \pkg{\link[BSgenome]{BSgenome-class}} object of reference genome to get reference bases from.
#' @importFrom methods is as
#' @importFrom Rsamtools indexBam scanBamFlag ScanBamParam FaFile scanFa scanFaIndex
#' @importFrom GenomicAlignments readGAlignments pileLettersAt
#' @importFrom GenomicRanges sort GRanges start end mcols
#' @importFrom GenomeInfoDb seqnames keepSeqlevels seqlevels
#' @importFrom BSgenome getSeq
#' @importFrom IRanges IRanges
#' @importFrom Biostrings DNAStringSet alphabetFrequency
#' @return A \code{list} containing two data frames objects called `mismatches` storing reported mismatches
#' and `coverage` reporting overall coverage per genomic position.
#' @author David Porubsky
#' @export
#' @examples
#' ## Get BAM to extract mismatches from
#' bam.file <- system.file("extdata", "test2.bam", package = "ggNucFreak")
#' ## Set region to extract mismatches from
#' region <- as('chr1:197373833-197391376', "GRanges")
#' ## Load bsgnome object to obtain reference bases from
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' ## Extract mismatches
#' m.obj <- exportMismatches(bamfile = bam.file, region = region, bsgenome = BSgenome.Hsapiens.UCSC.hg38)
#' m.obj
#'
exportMismatches <- function(bamfile=NULL, bamindex=bamfile, secondary.aln=FALSE, keep.duplicates=FALSE, region=NULL, min.mapq=10, min.baseq=20, min.mismatches=5, by.phase=FALSE, ref.fasta=NULL, bsgenome=NULL) {
  ## Helper function
  collapse.str <- function(x) {
    if (length(x) > 1) {
      paste(x, collapse='')
    } else if (length(x) == 1) {
      as.character(x)
    } else {
      x <- ""
    }
  }

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
  ## Get region positions
  base.pos <- GenomicRanges::GRanges(seqnames=GenomeInfoDb::seqnames(region.gr),
                                     ranges=IRanges::IRanges(start = seq( GenomicRanges::start(region.gr), GenomicRanges::end(region.gr), by=1 ),
                                                             end = seq( GenomicRanges::start(region.gr), GenomicRanges::end(region.gr), by=1) )
                                     )

  ## Get reference base at each position
  if (!is.null(bsgenome)) {
    ## Extract FASTA from BSgenome object
    ref.seq <- BSgenome::getSeq(bsgenome, region.gr)
  } else if (is.character(ref.fasta)) {
    ## Extract FASTA from user defined FASTA file
    fa.file <- open(Rsamtools::FaFile(ref.fasta))
    ## Remove sequences not present in the FASTA index
    fa.idx <- Rsamtools::scanFaIndex(fa.file)
    ## Read in sequence for a given range(s)
    ref.seq <- Rsamtools::scanFa(file = fa.file, param = region.gr, as = "DNAStringSet")
  } else {
    stop("Please set a 'bsgenome' or 'ref.fasta' parameter!!!")
  }
  ref.bases <- strsplit(as.character(ref.seq), '')[[1]]

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
  ## Filter by mapping quality
  if (min.mapq > 0) {
    aln.data <- aln.data[mcols(aln.data)$mapq >= min.mapq]
  }
  ## Sort alignments
  aln.data <- GenomicRanges::sort(aln.data, ignore.strand = TRUE)
  ## Restrict seqlevels to region of interest
  aln.data <- GenomeInfoDb::keepSeqlevels(aln.data, value = as.character(GenomeInfoDb::seqnames(region.gr)), pruning.mode = 'coarse')

  ## Get mismatches
  haps <- unique(GenomicRanges::mcols(aln.data)$HP)
  mism.l <- list()
  cov.l <- list()
  for (i in seq_along(haps)) {
    hap <- haps[i]
    sub.data <- aln.data[mcols(aln.data)$HP == hap]
    seqnames <- GenomeInfoDb::seqlevels(sub.data)
    ## Get nucleotide piles
    piles <- GenomicAlignments::pileLettersAt(sub.data@elementMetadata$seq, GenomeInfoDb::seqnames(sub.data), GenomicRanges::start(sub.data), sub.data@cigar, base.pos)
    ## Get total coverages
    coverage <- width(piles)
    ## Get base quality piles
    quals <- GenomicAlignments::pileLettersAt(sub.data@elementMetadata$qual, GenomeInfoDb::seqnames(sub.data), GenomicRanges::start(sub.data), sub.data@cigar, base.pos)
    ## Convert data frames
    df.piles <- methods::as(piles, "data.frame")
    df.quals <- methods::as(quals, "data.frame")
    ## Filter by user defined minimum base quality
    bases <- strsplit(df.piles[,1], "")
    quals <- sapply(df.quals[,1], function(x) as.numeric(charToRaw(x))-33)
    filtbases <- mapply(function(X,Y) { X[Y >= min.baseq] }, X=bases, Y=quals)
    ## Keep only non-reference bases per position
    nonref.bases <- mapply(function(X,Y) { X[X != Y] }, X=bases, Y=ref.bases)
    nonref.bases <- lapply(nonref.bases, collapse.str)
    ## Get non-reference allele frequency
    nonref.bases <- Biostrings::DNAStringSet(do.call(c, nonref.bases))
    nonref.bases.counts <- Biostrings::alphabetFrequency(nonref.bases)
    nonref.bases.counts <- nonref.bases.counts[, c('A', 'C', 'G', 'T')]
    ## Get mismatches
    #mism.count <- rowSums(nonref.bases.counts)
    #mask <- mism.count > 0
    mask <- apply(nonref.bases.counts, 1, function(x) any(x >= min.mismatches))
    ## Construct mismatch data.frame
    df <- as.data.frame(nonref.bases.counts[mask, , drop=FALSE])
    if (nrow(df) > 0) {
      df$chr <- seqnames
      df$pos <- GenomicRanges::start(base.pos)[mask]
      df <- reshape2::melt(df, measure.vars = c('A', 'C', 'G', 'T'))
      df$HP <- hap
    } else {
      df <- NULL
    }
    ## Get overall coverage per position
    df.cov <- data.frame(chr=rep(seqnames, length(ref.bases)), pos=start(base.pos), cov=coverage, HP=hap, ref.base=ref.bases)

    ## Store
    mism.l[[i]] <- df
    cov.l[[i]] <- df.cov
  }
  mism.df <- do.call(rbind, mism.l)
  cov.df <- do.call(rbind, cov.l)

  ## Return
  return(list(mismatches = mism.df, coverage = cov.df))
}

