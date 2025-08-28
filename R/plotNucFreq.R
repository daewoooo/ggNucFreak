#' Load BAM alignments a plot base frequencies.
#'
#' @param bamfile A path to a BAM file to get base frequencies from.
#' @param region A user defined region from where the base frequencies will be extracted.
#' @param max.data.points A maximum number of data points to be plotted for efficiency (Default : `100000`).
#' @param base.qual A minimum base quality to be 
#' @param min.mapq A minimum mapping quality of aligned reads to be extracted from provided BAM file.
#' @param max.depth (Default : `100`).
#' @param shape A user defined data point shape to visualize base frequency,either 'point' or 'step'.
#' @param nucleotide.content A user defined nucleotides whose frequency will be calculated, such 'GC' content.
#' @param bin.size A user defined 
#' @param reverse Set to \code{TRUE}
#' @param title A user defined title of the final plot.
#' @param return Define to return plot only `plot`, plotting data only `data` or both `c('plot', 'data')` (Default : `plot`).
#' @return A \code{ggplot} plot.
#' @import ggplot2
#' @importFrom methods as is
#' @importFrom scales comma
#' @importFrom deepSNV bam2R
#' @importFrom zoo rollapply
#' @importFrom S4Vectors Rle
#' @importFrom plyranges as_ranges
#' @importFrom GenomicRanges start end GRanges
#' @importFrom MatrixGenerics rowSums rowMaxs
#' @author David Porubsky
#' @export
#' @examples
#' ## Get BAM to extract methylation from
#' bam.file <- system.file("extdata", "test.bam", package = "ggNucFreak")
#' ## Set region to extract methylation from
#' region <- as('chr11:2787889-2790647', "GRanges")
#' region <- GenomicRanges::resize(region, width = GenomicRanges::width(region) + (GenomicRanges::width(region) * 2), fix = 'center')
#' ## Plot frequencies of maximum and the remainder base(s).
#' plotNucFreq(bamfile = bam.file, region = region)

plotNucFreq <- function(bamfile = NULL, region = NULL, max.data.points = 100000, base.qual = 25, min.mapq = 0, max.depth = 100, shape = 'point', nucleotide.content = NULL, bin.size = 100000, reverse = FALSE, title = NULL, return = 'plot') {
  ## Check user input
  if (!is.null(bamfile)) {
    if (!file.exists(bamfile)) {
      stop("BAM file defined in 'bamfile' does not exists !!!")
    }
  } else {
    stop("Parameter 'bamfile' point to the BAM file to process is required !!!")
  }

  if (!is.null(region)) {
    if (is.character(region)) {
      region.gr <- methods::as(region, "GRanges")
    } else if (methods::is(region, "GRanges")) {
      region.gr <- region
    }
  } else {
    stop("Parameter 'region' point to the genomic coordinates to process is required !!!")
  }

  ## Get nucleotide frequency ##
  region.size <- paste0(scales::comma(width(region.gr)), ' bp')
  ptm <- startTimedMessage("[deepSNV::bam2R] Reading in BAM alignments, region size: ", region.size, ' ...')
  pileup <- tryCatch(
    {
      deepSNV::bam2R(file = bamfile,
                     chr = as.character(seqnames(region.gr)),
                     start = as.numeric(start(region.gr)),
                     stop = as.numeric(end(region.gr)),
                     q = base.qual,
                     mq = min.mapq,
                     max.depth = max.depth)
    },
    error = function(e) {
      return(NULL)
    }  
  )
  stopTimedMessage(ptm)

  if (!is.null(pileup)) {
    ## Get nucleotide counts ##
    ###########################
    ptm <- startTimedMessage("[ggNucFreq] Constructing nucFreq plot ...")
    
    ## Reverse bases pileup if desired
    if (reverse) {
      pileup <- pileup[seq(dim(pileup)[1],1),]
    }
    
    plus.sum <- MatrixGenerics::rowSums(pileup[,c('A', 'T', 'C', 'G')])
    minus.sum <- MatrixGenerics::rowSums(pileup[,c('a', 't', 'c', 'g')])
    plus.max <- MatrixGenerics::rowMaxs(pileup[,c('A', 'T', 'C', 'G')])
    minus.max <- MatrixGenerics::rowMaxs(pileup[,c('a', 't', 'c', 'g')])
    plus.remainder <- plus.sum - plus.max
    minus.remainder <- minus.sum - minus.max
    max.total <- plus.max + minus.max
    remainder.total <- plus.remainder + minus.remainder
    total.bases <- max.total + remainder.total
    
    #plus.scnd <- apply(pileup[,c('A', 'T', 'C', 'G')], 1, function(x) sort(x,partial=4-1)[4-1])
    ## Get indels
    #indels.total <- rowSums(pileup[,c('INS', 'ins', 'DEL', 'del')])
  
    ## Get median coverage
    med.cov <- median(max.total + remainder.total)
  
    ## Reverse calculated vector if desired
    # if (reverse) {
    #   max.total <- rev(max.total)
    #   remainder.total <- rev(remainder.total)
    # }
  
    ## Get nucleotide content ##
    ############################
    if (!is.null(nucleotide.content) & methods::is(nucleotide.content, 'character')) {
      ## Split to separate letters and make sure all are valid bases
      letters <- strsplit(nucleotide.content, '')[[1]]
      if (all(grepl(letters, pattern='[A,C,G,T]', ignore.case = TRUE))) {
        letters2test <- c(toupper(letters), tolower(letters))
        dinuc.rowsums <- MatrixGenerics::rowSums(pileup[,letters2test])
        dinuc.bin.sums <- zoo::rollapply(data = dinuc.rowsums, width=bin.size, by=bin.size, FUN=sum)
        total.bin.sums <- zoo::rollapply(data = total.bases, width=bin.size, by=bin.size, FUN=sum)
        dinuc.bin.freq <- dinuc.bin.sums / total.bin.sums
        dinuc.bin.freq[is.nan(dinuc.bin.freq)] <- 0
        bin.starts <-  seq(from = as.numeric(GenomicRanges::start(region.gr)), to = (as.numeric(GenomicRanges::end(region.gr)) - bin.size) + 1, by = bin.size)
        bin.ends <- seq(from = (as.numeric(GenomicRanges::start(region.gr)) + bin.size) - 1, to = as.numeric(GenomicRanges::end(region.gr)), by = bin.size)
        dinuc.freq <- data.frame(start = bin.starts, end = bin.ends, bases = total.bin.sums, freq = dinuc.bin.freq)
        
        ## Plot nucleotide frequency ##
        ###############################
        ## Get color palette
        plt.data <- SVbyEye::getColorScheme(data.table = dinuc.freq, value.field = 'freq', breaks = c(0.25, 0.5, 0.75))
        ## Make a plot
        plt <- ggplot2::ggplot() +
          ggplot2::geom_rect(data=plt.data$data, ggplot2::aes(xmin=start, xmax=end, ymin=0, ymax=bases, fill=col.levels)) +
          ggplot2::scale_fill_manual(values = plt.data$colors, name=paste0(nucleotide.content, ' nucleotide\nfrequency'))  +
          ggplot2::scale_x_continuous(name = as.character(region.gr), expand = c(0, 0), label = comma) +
          ggplot2::scale_y_continuous(name = 'Total # of bases', label = scales::comma) +
          ggplot2::theme_minimal()
      }
    } else {
      ## Plot MAX and second max nucleotide frequency ##
      ##################################################
      ## Prepare base count table
      gen.pos <- seq(from = as.numeric(start(region.gr)), to = as.numeric(end(region.gr)))
      plt.data <- data.frame(gen.pos = gen.pos, max.base = max.total, remainder.base = remainder.total)
      ## Subset data to user defined max datapoints [default: 100000]
      filt <- sample(1:length(max.total))[1:max.data.points]
      max.df <- data.frame(gen.pos = gen.pos[filt], value = max.total[filt], Base = 'Max')
      max.df <- max.df[!is.na(max.df$value),]
      ## Compress remainder base counts by read-length-encoding
      remainder.rle <- S4Vectors::Rle(remainder.total)
      rle.ranges <- plyranges::as_ranges(remainder.rle)
      rle.pos <- c(rbind(start(rle.ranges), end(rle.ranges)))
      remainder.df <- data.frame(gen.pos = gen.pos[rle.pos],
                                 value = rep(runValue(remainder.rle), each=2),
                                 Base = 'Remainder'
      )
      remainder.df <- remainder.df[!is.na(remainder.df$value),]
  
      ## Make a plot
      plt.df <- rbind(max.df, remainder.df)
      if (shape == 'step') {
        plt <- ggplot2::ggplot() +
          ggplot2::geom_step(data = plt.df, ggplot2::aes(x=gen.pos, y = value, color = Base)) 
      } else if (shape == 'point') {
        plt <- ggplot2::ggplot() +
          ggplot2::geom_point(data = plt.df, ggplot2::aes(x=gen.pos, y = value, color = Base), size = 0.5)
      }
      plt <- plt +
        ggplot2::geom_hline(yintercept = med.cov, color='gold1', linetype='dashed', linewidth=1) +
        ggplot2::scale_color_manual(values = c('Max' = 'black', 'Remainder' = 'red'), name='Base count') +
        ggplot2::scale_x_continuous(expand = c(0, 0), label = scales::comma, name = as.character(region.gr)) +
        ggplot2::ylab('Max vs remainder base count') +
        ggplot2::theme_minimal()
    }
  
    ## Add plot title if defined
    if (!is.null(title)) {
      if (nchar( title) > 0) {
        plt <- plt + ggplot2::ggtitle(title)
      }
    }
  } else {
    plt <- ggplot2::ggplot()
    plt.data <- NULL
    message("Failed to extract read pileups from a BAM file ", bamfile, "\n", "corresponding to region ", as.character(region.gr), ", returning empty plot !!!")
  }  

  stopTimedMessage(ptm)
  ## Return final plot and/or data table
  if (all(return %in% c('plot', 'data'))) {
    return(list(plot = plt, data = plt.data))
  } else if (any(return %in% 'plot')) {
    return(plt) 
  } else if (any(return %in% 'data')) {
    return(plt.data) 
  } else {
    return(list(plot = plt, data = plt.data))
  }
}
