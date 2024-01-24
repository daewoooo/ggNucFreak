#' Collapse consecutive bins with the same ID value
#'
#' Collapse consecutive bins with the same value defined in 'id.field'. 
#'
#' @param gr A \code{\link{GRanges-class}} object.
#' @param id.field A number of metadata column to use for region merging.
#' @param measure.field A field column that contains measured value to be sum up.
#' @importFrom S4Vectors runLength Rle
#' @return A \code{\link{GRanges-class}} object.
#' @author David Porubsky
#' @export
#' 
collapseBins <- function(gr, id.field=0, measure.field=NULL) {
  ## Include seqnames into the unique ID
  unique.ID <- paste0(seqnames(gr), '_', mcols(gr)[,id.field])
  ## Get continous runs of the same unique ID
  unique.ID.runLength <- S4Vectors::runLength(S4Vectors::Rle(unique.ID))
  ind.last <- cumsum(unique.ID.runLength) ##get indices of last range in a consecutive(RLE) run of the same value
  ind.first <- c(1, cumsum(unique.ID.runLength) + 1) ##get indices of first range in a consecutive(RLE) run of the same value
  ind.first <- ind.first[-length(ind.first)]  ##erase last index from first range indices 
  collapsed.gr <- GenomicRanges::GRanges(seqnames=GenomeInfoDb::seqnames(gr[ind.first]), 
                                         ranges=IRanges::IRanges(start=start(gr[ind.first]), end=end(gr[ind.last])), 
                                         mcols=GenomicRanges::mcols(gr[ind.first]))
  names(mcols(collapsed.gr)) <- names(mcols(gr[ind.first]))
  ## Sum values in user defined measure fields(s)
  if (!is.null(measure.field)) {
    run.ID  <- factor(rep(seq_along(unique.ID.runLength), unique.ID.runLength))
    for (field in measure.field) {
      mcols(collapsed.gr)[,field] <-  as.numeric(tapply(mcols(gr)[,field], run.ID, sum))
      #collapsed.gr$C <-  tapply(gr$C, run.ID, sum)
      #collapsed.gr$W <-  tapply(gr$W, run.ID, sum)
    }
  }
  return(collapsed.gr)
}
