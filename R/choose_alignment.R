#' choose_alignment
#'
#' function to choose sample based on alignment using EMBL clustal omega
#'
#' @param lib_path string. path. location of strain library
#' @param lib_file string. file name. strain library file
#' @param blast_file string. file name. csv file of blast results
#' @param align_path string. path. location of alignment file
#' @param alignment string or vector of strings. file name. alignment file
#'     downloaded from EMBL clustal omega.
#' @param viewinR logical. default \code{FALSE.} When TRUE, displays alignment
#'     file in R console
#'
#' @return updates strain library with strains selected
#'
#' @import stringr
#' @export

choose_alignment <- function(lib_path, lib_file, blast_file,
                             align_path, alignment,
                             viewinR = FALSE) {
  # read in strain files-------------------------------------------------------
  lib_df <- read.csv(file.path(lib_path, lib_file), header = TRUE,
                     stringsAsFactors = FALSE)

  blast_df <- read.csv(blast_file, header = TRUE)

  for (i in 1:length(alignment)) {

    msg <- sprintf("\nReading alignment file: %s\n", alignment[i])
    message(msg)

    align_df <- read.table(file.path(align_path, alignment[i]),
                           sep = '\t', fill=TRUE, row.names = NULL)


    if(viewinR) {
      print(align_df)
    }

    # allowing user selection-----------------------------------------------------
    choices <- align_df[,1]
    choices <- unique(unlist(str_extract_all(choices, ".*ab1")))

    # identify the strain that is already in strain list
    if(any(choices %in% lib_df$sample_name)) {
      in_lib <- choices[choices %in% lib_df$sample_name]

      msg <- sprintf("\n%s is already found in strain library", in_lib)
      msg <- paste(msg, lib_file, sep = ' ')
      for(i in msg) message(msg)
    }


    prompt <- "Select which strains to add to strain library"
    selection <- select.list(title = prompt, choices = choices,
                             graphics = FALSE, multiple = TRUE)

    # adding blast result for selection to strain library
    entry <- blast_df[blast_df$sample_name %in% selection,]
    write.table(entry, file.path(lib_path, lib_file), append = TRUE, sep = ',', row.names = FALSE,
                col.names = FALSE)

    msg <- sprintf("Added blast result for %s to strain library",
                   selection)
    msg <- paste(msg, lib_file, sep = ' ')
    for(m in msg) message(m)
  }

  # making fasta file of strain library-----------------------------------------
  lib_name <- str_extract(lib_file, ".*(?=\\.csv)")
  msg <- sprintf("Writing fasta file of strain library: %s.fasta", lib_name)
  message(msg)

  lib_df <- read.csv(file.path(lib_path, lib_file), header = TRUE,
                     stringsAsFactors = FALSE)

  df2fasta(df =lib_df, fname = lib_name, path = lib_path,
           header_col = "sample_name", seq_col = "query_seq")
}
