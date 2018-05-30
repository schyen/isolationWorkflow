#' df2fasta
#'
#' converts dataframe into fasta file, where fasta headers are in one dataframe
#' column and sequences are in another dataframe column
#'
#' @param df dataframe
#' @param fname string. name of fasta file
#' @param path string. location of output fasta files
#' @param header_col column name containing sequence headers
#' @param seq_col column name of sequences


df2fasta <- function(df, fname = NULL, path, header_col, seq_col) {

  if(!header_col %in% colnames(df)) {
    stop("Value supplied to header_col is not found in dataframe")
  }

  if(!seq_col %in% colnames(df)) {
    stop("Value supplied to seq_col is not found in dataframe")
  }

  df <- df[,c(header_col, seq_col)]

  vert <- do.call(rbind, lapply(seq(nrow(df)), function(i) t(df[i, ])))

  out <- c()
  for(i in 1:nrow(vert)) {

    # header
    if(i %% 2 == 1) {
      entry <- sprintf(">%s", vert[i])
      entry <- str_replace_all(entry, ' ', '-')
    }
    else {
      entry <- vert[i]
    }

    out <- rbind(out, entry)
  }

  write.table(out, file = file.path(path, sprintf("%s.fasta", fname)),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}
