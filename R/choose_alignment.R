#' choose_alignment
#'
#' function to choose sample based on alignment using EMBL clustal omega
#'
#' @param lib_path string. path. location of strain library
#' @param lib_file string. file name. strain library file
#' @param blast_file string. file name. csv file of blast results
#' @param fasta_path string. path. location of fasta file(s) that need to be aligned
#' @param fastatoalign string or vector of strings. fasta files to be aligned.
#'
#' @return updates strain library with strains selected
#'
#' @import stringr
#' @export

choose_alignment <- function(lib_path, lib_file, blast_file,
                             fasta_path, fastatoalign) {
  # read in strain files-------------------------------------------------------
  lib_df <- read.csv(file.path(lib_path, lib_file), header = TRUE,
                     stringsAsFactors = FALSE)

  blast_df <- read.csv(blast_file, header = TRUE)

  for (i in 1:length(fastatoalign)) {
    print(fastatoalign[i])
    # reading in sequences
    str_set <- Biostrings::readDNAStringSet(file.path(fasta_path, fastatoalign[i]))

    # remove gaps
    str_set <- DECIPHER::RemoveGaps(str_set, removeGaps = "all")

    # perform alignment
    align_file <- str_extract(fastatoalign[i], ".*(?=\\.fasta)")
    align_file <- paste0("alignment-", align_file, ".html")
    align_file <- file.path(fasta_path, align_file)

    print(align_file)
    aligned <- DECIPHER::AlignSeqs(str_set, iterations = 0, refinements = 0)

    # view alignment
    DECIPHER::BrowseSeqs(aligned, htmlFile = align_file)

    # sample names
    sample_name <- names(aligned)

    # allowing user selection-----------------------------------------------------
    # choices <- align_df[,1]
    choices <- sample_name
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
