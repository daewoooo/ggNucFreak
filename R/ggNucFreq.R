## Function definitions ##
##########################

#' Load BAM alignments a plot base frequencies.
#'
#' This function takes loaded PAF alignments using \code{\link{readPaf}} function and perform
#' user defined filtering of input alignments based on mapping quality, alignment length, and
#' minimum alignments between target and query.
#'
#' @param bam.file ...
#' @param region ...
#' @param max.data.points ...
#' @param base.qual ...
#' @param map.qual ...
#' @param max.depth ...
#' @param shape ...
#' @param nucleotide.content ...
#' @param bin.size ...
#' @param reverse Set to \code{TRUE}
#' @param title ...
#' @return A \code{ggplot} plots
#' @importFrom dplyr group_by mutate n
#' @importFrom utils read.table
#' @author David Porubsky
#' @export
ggNucFreq <- function(bam.file = NULL, region = NULL, max.data.points = 100000, base.qual = 25, map.qual = 0, max.depth = 100, shape = 'step', nucleotide.content = NULL, bin.size = 100000, reverse = FALSE, title = NULL) {
  ## Check user input
  if (!is.null(bam.file)) {
    if (!file.exists(bam.file)) {
      stop("BAM file defined in 'bam.file' does not exists !!!")
    }
  } else {
    stop("Parameter 'bam.file' point to the BAM file to process is required !!!")
  }

  if (!is.null(region)) {
    if (is.character(region)) {
      region.gr <- methods::as(region, "GRanges")
    } else if (is(region, "GRanges")) {
      region.gr <- region
    }
  } else {
    stop("Parameter 'region' point to the genomic coordinates to process is required !!!")
  }

  ## Get nucleotide frequency ##
  region.size <- paste0(scales::comma(width(region.gr)), ' bp')
  ptm <- startTimedMessage("[deepSNV::bam2R] Reading in BAM alignments, region size: ", region.size, ' ...')
  pileup <- deepSNV::bam2R(file = bam.file,
                           chr = as.character(seqnames(region.gr)),
                           start = as.numeric(start(region.gr)),
                           stop = as.numeric(end(region.gr)),
                           q = base.qual,
                           mq = map.qual,
                           max.depth = max.depth
  )
  stopTimedMessage(ptm)

  ## Get nucleotide counts ##
  ###########################
  ptm <- startTimedMessage("[ggNucFreq] Constructing nucFreq plot ...")
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
  if (reverse) {
    max.total <- rev(max.total)
    remainder.total <- rev(remainder.total)
  }

  ## Get dinucleotide content ##
  #############################
  if (!is.null(nucleotide.content) & methods::is(nucleotide.content, 'character')) {
    ## Split to separate letter and make sure all are valid bases
    letters <- strsplit(nucleotide.content, '')[[1]]
    if (all(grepl(letters, pattern='[A,C,G,T]', ignore.case = TRUE))) {
      letters2test <- c(toupper(letters), tolower(letters))
      dinuc.rowsums <- MatrixGenerics::rowSums(pileup[,letters2test])
      dinuc.bin.sums <- zoo::rollapply(data = dinuc.rowsums, width=bin.size, by=bin.size, FUN=sum)
      total.bin.sums <- zoo::rollapply(data = total.bases, width=bin.size, by=bin.size, FUN=sum)
      dinuc.bin.freq <- dinuc.bin.sums / total.bin.sums
      dinuc.bin.freq[is.nan(dinuc.bin.freq)] <- 0
      bin.starts <-  seq(from = as.numeric(start(region.gr)), to = (as.numeric(end(region.gr)) - bin.size) + 1, by = bin.size)
      bin.ends <- seq(from = (as.numeric(start(region.gr)) + bin.size) - 1, to = as.numeric(end(region.gr)), by = bin.size)
      dinuc.freq <- data.frame(start = bin.starts, end = bin.ends, bases = total.bin.sums, freq = dinuc.bin.freq)

      ## Plot dinucleotide frequency ##
      #################################
      ## Get color palette
      plt.data <- SVbyEye::getColorScheme(data.table = dinuc.freq, value.field = 'freq', breaks = c(0.25, 0.5, 0.75))
      ## Make a plot
      plt <- ggplot() +
        geom_rect(data=plt.data$data, aes(xmin=start, xmax=end, ymin=0, ymax=bases, fill=col.levels)) +
        scale_fill_manual(values = plt.data$colors, name=paste0(nucleotide.content, ' nucleotide\nfrequency'))  +
        scale_x_continuous(name = 'Genomic position (bp)', expand = c(0, 0), label = comma) +
        scale_y_continuous(name = 'Total # of bases', label = comma) +
        theme_minimal()
    }
  } else {
    ## Plot MAX and second max nucleotide frequency ##
    ##################################################
    ## Subset data to user defined max datapoints [default: 100000]
    gen.pos <- seq(from = as.numeric(start(region.gr)), to = as.numeric(end(region.gr)))
    filt <- sample(1:length(max.total))[1:max.data.points]
    max.df <- data.frame(gen.pos = gen.pos[filt], max.base = max.total[filt])
    ## Compress remainder base counts by read-length-encoding
    remainder.rle <- S4Vectors::Rle(remainder.total)
    rle.ranges <- plyranges::as_ranges(remainder.rle)
    rle.pos <- c(rbind(start(rle.ranges), end(rle.ranges)))
    remainder.df <- data.frame(gen.pos = gen.pos[rle.pos],
                               value = rep(runValue(remainder.rle), each=2)
    )

    ## Make a plot
    if (shape == 'step') {
      plt <- ggplot() +
        geom_step(data = max.df, aes(x=gen.pos, y = max.base, color = 'Max.base')) +
        geom_step(data = remainder.df, aes(x=gen.pos, y = value, color = 'Remainder.base(s)'))
    } else if (shape == 'point') {
      plt <- ggplot() +
        geom_point(data = max.df, aes(x=gen.pos, y = max.base, color = 'Max.base'), size = 0.5) +
        geom_point(data = remainder.df, aes(x=gen.pos, y = value, color = 'Remainder.base(s)'), size = 0.5)
    }
    plt <- plt +
      geom_hline(yintercept = med.cov, color='gold1', linetype='dashed', linewidth=1) +
      scale_color_manual(values = c('Max.base' = 'black', 'Remainder.base(s)' = 'red'), name='') +
      scale_x_continuous(expand = c(0, 0), label = scales::comma, name = 'Genomic position (bp)') +
      ylab('Max vs remainder base count') +
      theme_minimal()
  }

  ## Add plot title if defined
  if (!is.null(title)) {
    if (nchar( title) > 0) {
      plt <- plt + ggtitle(title)
    }
  }

  ## Return results
  stopTimedMessage(ptm)
  return(plt)
}



# ## Load required libraries
# library(deepSNV)
# library(ggplot2)
# library(Rsamtools)
#
# ## Load region index file
# infile <- '/home/porubsky/WORK/Projects/ggNucFreq/chr22_18000000_23000000_index.tsv'
# idx.df <- read.table(infile, header = TRUE, sep = '\t', stringsAsFactors = FALSE, comment.char = '&')
# ## Keep only HPRC samples
# idx.df <- idx.df[grep(idx.df$fasta.file, pattern='pat|mat'),]
#
# plots <- list()
# for (i in 1:nrow(idx.df)) {
#   loc.df <- idx.df[i,]
#   region <- makeGRangesFromDataFrame(loc.df)
#   ## Get sample id
#   hap.id <- loc.df$hap.id
#   sample.id <- gsub(hap.id, pattern = '_1|_2', replacement = '')
#   ## Get corresponding alignments
#   bam.path <- file.path('/home/porubsky/EEE_volumes/vol27/projects/hprc/nobackups/analysis/sedef/realign/', sample.id, '/alignments/HiFi/minimap2/')
#   bam.file <- list.files(bam.path, pattern = '\\.bam$', full.names = TRUE)
#
#   plt <- ggNucFreq(bam.file = bam.file, region = region, reverse = loc.df$revcomp,  title = hap.id)
#   plots[[i]] <- plt
# }
#
# ## Testing ###################################################################################
# plt.df <- data.frame(gen.pos = gen.pos, max.base = max.total, remainder = remainder.total)
# plt.df <- plt.df[sample(1:nrow(plt.df))[1:100000],]
#
# remainder.rle <- Rle(plt.df$remainder)
# rle.ranges <- as_ranges(remainder.rle)
# rle.pos <- c(rbind(start(rle.ranges), end(rle.ranges)))
# remainder.df <- data.frame(gen.pos = plt.df$gen.pos[rle.pos],
#                            value = rep(runValue(remainder.rle), each=2)
# )
#
# ggplot() +
#   geom_step(data = plt.df, aes(x=gen.pos, y = max.base), color = 'black') +
#   #geom_point(data =  plt.df, aes(x=gen.pos, y = remainder), color = 'red') +
#   geom_step(data = remainder.df, aes(x=gen.pos, y = value), color = 'red') +
#   theme_minimal()
#
#
# ptm <- startTimedMessage("Rsamtools")
# sbp <- ScanBamParam(which=GRanges(loc.df$seqnames, IRanges(loc.df$start, loc.df$end)))
# p_param <- PileupParam(min_base_quality = 10L)
# res <- pileup(bam.file, scanBamParam=sbp, pileupParam = p_param)
# stopTimedMessage(ptm)
