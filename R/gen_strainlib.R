#' gen_strainlib
#'
#' function to generate/update strain library and fasta file
#' reads in blast result file from parse blast
#' if no strain library supplied, creates a new one
#' Species that are only found once and not already in strain library are
#' added to library
#' Generates strain library as csv
#' Generates strain library as fasta
#' also generates fasta file for each species, amalgamating results from
#' current blast result and strain library
#'
#' @param blast_file string. file name. csv file containing blast results. Include path if necessary
#' @param lib_path path. string. location of strain library file
#' @param lib_file string. file name. defualt NULL. Name of strain library file.
#'     Creates a library file if none supplied.
#' @param lib_name string. defaults to generic library name,name of strain library (only used if new one is
#'     being created)
#'
#' @return creates/updates strain library
#'     creates fasta file of strain library
#'     creates fasta files of species with multiple isolates found
#'
#' @import dyplr
#' @import stringr
#' @export

gen_strainlib <- function(blast_file, lib_path, lib_file = NULL,
                               lib_name = 'mystrainlibrary') {

  # read in blast result--------------------------------------------------------
  blast <- read.csv(blast_file, header = TRUE)

  # if strain library not supplied, make new, empty library---------------------
  if(is.null(lib_file)) {

    lib_file <- paste0(file.path(lib_path, lib_name), ".csv")

    # empty library df
    empty_lib <- blast[0, ]

    write.table(empty_lib, lib_file, col.names = TRUE, row.names = FALSE,
                sep = ',', na = "")
    msg <- sprintf("No strain library supplied. Creating new strain library: \n%s\n",
                   lib_file)
    message(msg)
  }
  else {
    lib_file <- file.path(lib_path, lib_file)
  }

  # read in library------------------------------------------------------------
  lib_df <- read.csv(lib_file, header = TRUE, stringsAsFactors = FALSE)


  # identify single species-----------------------------------------------------
  species_tally <- as.data.frame(table(blast$match))
  colnames(species_tally) <- c('match','Freq')

  single <- species_tally[species_tally$Freq == 1,]

  add_singles <- !single %in% lib_df$match

  add_singles <- single[add_singles]$match

  msg <- sprintf("%s found as unique species that is not currently in your strain library. Adding to strain library.", add_singles)
  for(i in msg) message(i)

  entry <- blast[blast$match == add_singles,]

  # appending entry to library
  write.table(entry, lib_file, append = TRUE, sep = ',', row.names = FALSE,
              col.names = FALSE)

  # Now go through rest of sequences and create fasta files for each species----
  # getting entries from blast result
  wip_df <- blast[!blast$match %in% add_singles,]

  # getting entries from strain library
  wip_df <- rbind(wip_df, lib_df[!lib_df$match %in% add_singles,])

  msg <- sprintf("\nGoing through current blast file and strain library, multiple strains of %s were found. Generating fasta file for this species for sequence alignment", unique(wip_df$match))

  for(i in msg) message(i)

  species_df <- wip_df[,c('sample_name','match','query_seq')] %>%
    group_by(match) %>%
    do(df2fasta(., fname = unique(.$match), path = lib_path,
                header_col = 'sample_name', seq_col = 'query_seq'))


}
