#' Retrieve Genomic Coordinates for a Specific Gene
#'
#' This function returns the genomic coordinates associated with a specified Ensembl gene ID.
#' The function uses GRCh37 genome annotations to obtain this information.
#' It fetches the coordinates of the most highly expressed transcript of the gene,
#' including its exons and coding sequence (CDS), and calculates the total length
#' of the transcript. The function filters out genes that are not located on one of the
#' standard chromosomes (1-22, X, Y).
#'
#' @importFrom data.table data.table
#' @param ensembl_gene_id An Ensembl gene identifier
#'
#' @return A data table object with coordinates
#' @export
#'
#' @examples
#' \dontrun{
#' ensembl_gene_id <- "ENSG00000133703"
#' exonsCoordinates <- T
#' get_gene_location(ensembl_gene_id, exonsCoordinates = F)
#' }
get_gene_location <- function(ensembl_gene_id) {
  # Create objects for exon coordinates and genes coordinates from annotation for GRCh37
  gtf_exons <- mutmatch_gtf_grch37[mutmatch_gtf_grch37$type == "exon"]
  gtf_genes <- mutmatch_gtf_grch37[mutmatch_gtf_grch37$type == "gene"]
  # gtf_cds <- mutmatch_gtf_grch37[mutmatch_gtf_grch37$type == "CDS"]
  gtf_stop <- mutmatch_gtf_grch37[mutmatch_gtf_grch37$type == "stop_codon"]
  # "exons" don't contain UTRs in my case since I filter  by CDS
  gtf_cds <- mutmatch_gtf_grch37[mutmatch_gtf_grch37$type %in% c("CDS", "stop_codon")]

  # Get chromosome coordinates for a gene (includes UTRs)
  gene_coordinates <- gtf_genes[gtf_genes$gene_id %in% ensembl_gene_id]
  gene_coordinates <- GenomicRanges::union(gene_coordinates, gene_coordinates)
  gene_coordinates <- data.table::as.data.table(gene_coordinates)[, -c(4)] %>% data.table::data.table()
  base::colnames(gene_coordinates) <- c("chromosome", "gene_start", "gene_end", "strand")

  # Only get genes on chromosomes
  gene_coordinates <- gene_coordinates[chromosome %in% c(1:22, "X", "Y")]

  # Retrieve coordinates for top-5 most expressed transcripts for a gene
  corresponding_transcripts <- mutmatch_topTranscripts %>%
    dplyr::filter(Gene %in% ensembl_gene_id) %>%
    dplyr::pull(Transcript_ID) %>%
    base::unique()
  transcript_info <- gtf_cds[gtf_cds$transcript_id %in% corresponding_transcripts]

  # If GTF CDS annotation does not have any transcript
  # from top-5 most expressed, take exons
  if (length(transcript_info) == 0) {
    transcript_info <- gtf_exons[gtf_exons$transcript_id %in% corresponding_transcripts]
  }

  # Obtain the most expressed transcript for a gene
  transcript_info <- base::merge(transcript_info, mutmatch_topTranscripts[, c("Transcript_ID", "Rank")],
    by.x = "transcript_id",
    by.y = "Transcript_ID"
  )
  transcript_info <- transcript_info %>%
    dplyr::filter(Rank == base::min(Rank)) %>%
    dplyr::select(-Rank)

  transcript_info <- GenomicRanges::makeGRangesFromDataFrame(transcript_info, keep.extra.columns = T)
  transcript_info <- GenomicRanges::union(transcript_info, transcript_info)
  transcript_info <- data.table::as.data.table(transcript_info)[, c(1:3)]
  colnames(transcript_info) <- c("chromosome", "exon_start", "exon_end")

  # Write the length of the gene transcript
  transcript_length <- base::sum((transcript_info$exon_end - transcript_info$exon_start) + 1,
    na.rm = T
  )
  gene_coordinates$length <- transcript_length

  gene_coordinates <- stats::na.omit(base::merge(gene_coordinates, transcript_info))
  return(gene_coordinates)
}
