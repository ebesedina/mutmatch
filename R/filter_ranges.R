#' Filter Genomic Ranges#'
#'
#' This function performs comprehensive filtering of genomic ranges, always excluding regions with low mappability
#' and conversion-unstable positions.
#' Additionally, it can optionally incorporate user-defined regions for customized filtering.
#'
#' @param grObject GenomicRanges object to be filtered.
#' @param filterRegionsPath File path to the additional regions that need to be filtered. Default is NULL.
#' @param format Format of the filterRegions file, can be "Rdata", "bed", or "NA".
#' "NA" means no additional filtering is performed. Default is "NA".
#' @param filter Indicates whether to include ("in") or exclude ("out") the regions specified in filterRegionsPath.
#' "NA" means no additional filtering is performed. Default is "NA".
#'
#' @return A GenomicRanges object after applying the specified filters.
#' @export
#' @examples
#' \dontrun{
#' filterRegionsPath <- system.file("extdata",
#'   "HotSpotAnnotations_Trevino.RData",
#'   package = "mutmatch"
#' )
#' }
filter_ranges <- function(grObject,
                          filterRegionsPath = NULL,
                          format = c("NA", "Rdata", "bed"),
                          filter = c("NA", "in", "out")) {
  # Validate arguments
  if (is.na(format)) {
    format <- "NA"
  } else {
    format <- base::match.arg(format)
  }

  if (is.na(filter)) {
    filter <- "NA"
  } else {
    filter <- base::match.arg(filter)
  }

  # Import and apply low mappability filter
  CRG_mappability_path <- system.file("extdata", "wgEncodeCrgMapabilityAlign75merMask.bed.gz", package = "mutmatch")
  CRG_mappability <- BiocIO::import(con = CRG_mappability_path, which = grObject, format = "bed")
  grObject <- GenomicRanges::intersect(grObject, CRG_mappability, ignore.strand = TRUE)

  # Import and apply conversion-unstable positions filter
  CUP_path <- system.file("extdata", "FASTA_BED.ALL_GRCh37.novel_CUPs.bed.gz", package = "mutmatch")
  CUP <- BiocIO::import(con = CUP_path, which = grObject, format = "bed")
  grObject <- GenomicRanges::setdiff(grObject, CUP, ignore.strand = TRUE)

  # Initialize variable for additional regions
  regions_mask <- NULL

  # Check for additional regions and import based on format
  if (!is.null(filterRegionsPath)) {
    if (format == "Rdata") {
      regions_mask <- base::readRDS(file = filterRegionsPath)
    } else if (format == "bed") {
      regions_mask <- BiocIO::import(con = filterRegionsPath, format = "bed")
    } else if (format == "NA") {
      message("No additional filtering is applied as format is set to 'NA'.")
    } else {
      stop("The 'format' parameter should be either 'Rdata', 'bed', or 'NA'+")
    }
  }

  # Apply additional filtering based on 'filter' argument
  if (!is.null(regions_mask) && filter != "NA") {
    if (filter == "in") {
      grObject <- GenomicRanges::intersect(grObject, regions_mask, ignore.strand = TRUE)
    } else if (filter == "out") {
      grObject <- GenomicRanges::setdiff(grObject, regions_mask, ignore.strand = TRUE)
    } else {
      stop("The 'filter' parameter should be either 'in' or 'out' if 'filterRegionsPath' is not NULL.")
    }
  }

  # Union to merge overlapping or duplicate ranges
  grObject <- GenomicRanges::reduce(grObject)

  return(grObject)
}
