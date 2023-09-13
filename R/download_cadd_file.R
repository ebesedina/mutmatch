#' Download Large CADD File with Real-time Progress Indicator
#'
#' Downloads a substantial CADD (Combined Annotation-Dependent Depletion) file
#' from its official repository (https://krishna.gs.washington.edu/download/CADD/bigWig/).
#' The file is saved to a designated location on the user's local file system.
#' Utilizes the Unix `wget` command-line utility to manage the download process and
#' to provide a real-time progress update.
#'
#' **Warning**: The target file is sizeable (~11GB); ensure adequate storage space and
#' anticipate extended download time. Windows users might experience problems with this command. A possible solution
#' would be to directly download the file from https://krishna.gs.washington.edu/download/CADD/bigWig/CADD_GRCh37-v1.4.bw.
#'
#' @param caddScoresPath A character string representing the full path where the
#' downloaded file should be saved, including the filename.
#'
#' @return NULL. The function executes for its side effects, which include downloading
#' the file and saving it to `caddScoresPath`.
#'
#' @examples
#' \dontrun{
#' # Replace "your/destination/path/CADD_GRCh37-v1.4.bw" with the path
#' # where you want to save the file.
#'
#' download_cadd_file(caddScoresPath = "your/destination/path/CADD_GRCh37-v1.4.bw")
#' }
#'
#' @export
download_cadd_file <- function(caddScoresPath) {
  # Specifying the URL where the CADD file can be downloaded.
  url <- "https://krishna.gs.washington.edu/download/CADD/bigWig/CADD_GRCh37-v1.4.bw"

  # Use `system` to invoke the wget command for downloading the file.
  # The `-O` option specifies the location where the file will be saved.
  system(paste("wget", url, "-O", caddScoresPath))

  # Displaying a message to indicate where the file has been downloaded and saved.
  cat("File successfully downloaded to:", caddScoresPath, "\n")
}
