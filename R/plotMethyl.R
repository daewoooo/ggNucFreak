#' Function to plot methylated cytosines from the BAM file in forward reference coordinates.
#'
#' This function takes loaded methylation porfiles using \code{\link{extractMethylCs}} function and visualize them 
#' using basic \pkg{ggplot2} functions. 
#'
#' @param by.phase Set to \code{TRUE} if the methylation should be plotted separated by haplotypes.
#' @param highlight.pos A vector of x-axis coordinates to be highlighted by the a solid line.
#' @param highlight.gr A \code{\link{GRanges-class}} object of all regions to be highlighted by a transparent rectangle.
#' @inheritParams extractMethylCs
#' @import ggplot2
#' @importFrom GenomicRanges start end
#' @importFrom IRanges subsetByOverlaps
#' @return A \code{\link{GRanges-class}} object with mirrored coordinates.
#' @author David Porubsky
#' @export
#' @examples
#' ## Get BAM file to plot
#' bam.file <- system.file("extdata", "test.bam", package = "ggNucFreak")
#' ## Set region to plot methylation from
#' region <- as('chr11:2787889-2790647', "GRanges")
#' region.plus <- GenomicRanges::resize(region, width = GenomicRanges::width(region) + (GenomicRanges::width(region) * 2), fix = 'center')
#' ## Plot methylation
#' plotMethyl(bamfile = bam.file, region = region.plus)
#' ## Plot methylation by phase
#' plotMethyl(bamfile = bam.file, region = region.plus, by.phase = TRUE)
#' ## Highlight specific positions
#' plotMethyl(bamfile = bam.file, region = region.plus, by.phase = TRUE, highlight.pos = c(2787889, 2790647))
#' ## Highlight specific region
#' plotMethyl(bamfile = bam.file, region = region.plus, by.phase = TRUE, highlight.gr = region)
#'
plotMethyl <- function(bamfile=NULL, region=NULL, min.mapq=10, by.phase=FALSE, highlight.pos=NULL, highlight.gr=NULL) {
  ## Check user input
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

  ## Get methylation data
  if (by.phase) {
    mCs <- extractMethylCs(bamfile = bamfile, region=region.gr, add.phase = TRUE)
  } else {
    mCs <- extractMethylCs(bamfile = bamfile, region=region.gr)
  }
  methylCs.gr <- mCs$methylCs
  reads.gr <- mCs$reads

  ## Set plotting theme
  custom_theme <- ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                                 axis.text.y = ggplot2::element_blank(),
                                 axis.ticks.y = ggplot2::element_blank(),
                                 axis.line.x = ggplot2::element_line(colour = "black", linewidth = 1, linetype = "solid"),
                                 axis.ticks.x = ggplot2::element_line(colour = "black", linewidth = 1, linetype = "solid"),
                                 panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank())

  ## Keep only C's with reported methylation
  plt.gr <- methylCs.gr[!is.na(methylCs.gr$m.prob)]

  ## Prepare data for plotting
  plt.df <- as.data.frame(plt.gr)
  ## Mask certain columns
  mask.columns <- which(colnames(GenomicRanges::mcols(reads.gr)) %in% c('seq', 'MM', 'ML'))
  reads.df <- as.data.frame(reads.gr[,-mask.columns])

  ## Make a plot
  if (by.phase) {
    reads.df$HP <- paste('Haplotype', reads.df$HP, sep = ' ')
    plt.df$HP <- paste('Haplotype', plt.df$HP, sep = ' ')
    plt <- ggplot2::ggplot() +
      ggplot2::geom_rect(data = reads.df, ggplot2::aes(xmin=start.crop, xmax=end.crop, ymin=phase.level - 1, ymax=phase.level, fill=strand)) +
      ggplot2::geom_linerange(data = plt.df, ggplot2::aes(x=start, ymin=phase.level - 1, ymax=phase.level, color=m.prob)) +
      ggplot2::geom_rect(data = reads.df, ggplot2::aes(xmin=start.crop, xmax=end.crop, ymin=phase.level - 1, ymax=phase.level), color='black', fill=NA) +
      ggplot2::facet_grid(HP ~ ., space = 'free', scales = 'free')
  } else {
    plt <- ggplot2::ggplot() +
      ggplot2::geom_rect(data = reads.df, ggplot2::aes(xmin=start.crop, xmax=end.crop, ymin=level - 1, ymax=level, fill=strand)) +
      ggplot2::geom_linerange(data = plt.df, ggplot2::aes(x=start, ymin=level - 1, ymax=level, color=m.prob)) +
      ggplot2::geom_rect(data = reads.df, ggplot2::aes(xmin=start.crop, xmax=end.crop, ymin=level - 1, ymax=level), color='black', fill=NA)
  }

  ## Add color scheme and plotting theme
  plt <- plt +
    ggplot2::scale_x_continuous(name='Genomic position (bp)') +
    ggplot2::scale_color_gradient(low = 'blue', high = 'red') +
    ggplot2::scale_fill_manual(values = c('+' = 'gray90', '-' = 'gray80')) +
    ggplot2::theme_minimal() + custom_theme

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

