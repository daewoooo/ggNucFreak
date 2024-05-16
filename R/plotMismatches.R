#' Function to plot positions of mismatched bases extracted from the BAM file using \code{\link{exportMismatches}} function.
#'
#' @param mismatches.obj \code{list} containing reported `mismatches` and `coverage` per genomic position.
#' @param add.reference.base Set to \code{TRUE} if reference base should be added to the plot.
#' @param coverage.profile Visulize overall coverage either as a solid 'line' or filled area 'fill'.
#' @inheritParams plotMethyl
#' @import ggplot2
#' @importFrom methods is as
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges subsetByOverlaps
#' @importFrom GenomeInfoDb seqnames keepSeqlevels seqlevels
#' @importFrom BSgenome getSeq
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
#' ## Extract and plot mismatches
#' m.obj <- exportMismatches(bamfile = bam.file, region = region, bsgenome = BSgenome.Hsapiens.UCSC.hg38)
#' plotMismatches(mismatches.obj = m.obj)
#' ## Extract and plot mismatches per haplotype
#' m.obj <- exportMismatches(bamfile = bam.file, region = region, bsgenome = BSgenome.Hsapiens.UCSC.hg38, by.phase=TRUE)
#' plotMismatches(mismatches.obj = m.obj)
#'
plotMismatches <- function(mismatches.obj, add.reference.base=FALSE, coverage.profile='line', highlight.pos=NULL, highlight.gr=NULL) {
  ## Check user input ##
  if (!methods::is(mismatches.obj, 'list')) {
    stop("Submitted object has to be a list containing data table 'mismatches' and 'coverage' exported by 'exportMismatches'function !!!")
    if (!all(names(mismatches.obj) %in% c("mismatches", "coverage"))) {
      stop("Submitted object does not contain expected data tables 'mismatches' and 'coverage' exported by 'exportMismatches'function !!!")
    }
  }

  ## Get data to plot
  mism.df <- mismatches.obj$mismatches
  cov.df <- mismatches.obj$coverage

  ## Add haplotype labels
  if (!is.null(mism.df)) {
    mism.df$HP <- ifelse(mism.df$HP == 0, 'Unphased', paste0('Haplotype ', mism.df$HP))
  }
  if (!is.null(cov.df)) {
    cov.df$HP <- ifelse(cov.df$HP == 0, 'Unphased', paste0('Haplotype ', cov.df$HP))
  }

  ## Get covered region
  cov.gr <- GenomicRanges::GRanges(seqnames = cov.df$chr, ranges=IRanges::IRanges(start=cov.df$pos, end = cov.df$pos), cov=cov.df$cov, HP=cov.df$HP)
  region.gr <- range(cov.gr)

  ## Make a plot ##
  #################
  ## Plot coverage profile
  if (coverage.profile == 'line') {
    plt <- ggplot2::ggplot() +
      ggplot2::geom_step(data = cov.df, ggplot2::aes(x=pos - 0.5, y=cov))
  } else if (coverage.profile == 'fill') {
    #cov.gr <- collapseBins(sort(cov.gr), id.field = 1)
    cov.ranges.df <- as.data.frame(cov.gr)

    plt <- ggplot2::ggplot() +
      ggplot2::geom_rect(data = cov.ranges.df, ggplot2::aes(xmin=start, xmax=end, ymin=0, ymax=cov), color='gray90', fill='gray90')
  }

  ## Add mismatches
  if (!is.null(mism.df)) {
    plt <- plt +
      ggplot2::geom_col(data = mism.df, ggplot2::aes(x=pos, y=value, fill=variable, color=variable), position = ggplot2::position_stack(), width = 0.75)
  }

  ## Add scales and faceting
  plt <- plt +
    ggplot2::scale_color_manual(values = c('A' = 'coral2', 'C' = 'lightgreen', 'G' = 'steelblue3', 'T' = 'tan2'), name='Base') +
    ggplot2::scale_fill_manual(values = c('A' = 'coral2', 'C' = 'lightgreen', 'G' = 'steelblue3', 'T' = 'tan2'), name='Base') +
    ggplot2::scale_x_continuous(name='Genomic position (bp)') +
    ggplot2::scale_y_continuous(name='Coverage') +
    ggplot2::facet_grid(HP ~ ., space='free', scale='free') +
    ggplot2::theme_minimal()

  ## Add reference base if desired
  if (add.reference.base) {
    ## Get approximate track width
    track.width <- max(cov.df$cov) * 0.1
    plt <- plt +
      ggplot2::geom_rect(data = cov.df, ggplot2::aes(ymin=0, ymax=-track.width, xmin=pos - 0.5, xmax=pos + 0.5, color=ref.base, fill=ref.base))
  }

  ## Highlight positions
  if (!is.null(highlight.pos)) {
    highlight.pos <- highlight.pos[highlight.pos >= GenomicRanges::start(region.gr) & highlight.pos <= GenomicRanges::end(region.gr)]
    plt <- plt + ggplot2::geom_vline(xintercept = highlight.pos)
  }

  ## Highlight regions
  if (!is.null(highlight.gr)) {
    if (methods::is(highlight.gr, 'GRanges')) {
      highlight.gr <- IRanges::subsetByOverlaps(highlight.gr, region)
      highlight.df <- as.data.frame(highlight.gr)
      plt <- plt + ggplot2::geom_rect(data=highlight.df, ggplot2::aes(xmin=start, xmax=end, ymin=0, ymax=Inf), fill='burlywood2', alpha=0.3)
    }
  }

  ## Return final plot
  return(plt)
}
