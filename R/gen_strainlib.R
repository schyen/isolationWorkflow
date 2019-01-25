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
#' @import dplyr
#' @import stringr
#' @export

gen_strainlib <- function(blast_file, lib_path, lib_file = NULL,
                               lib_name = 'mystrainlibrary') {

  # Check input files
  if(!file.exists(blast_file)) {
    stop(sprintf("\nFile does not found: %s", blast_file))
  }
  if(!dir.exists(lib_path)) {
    stop(sprintf("\nLocation not found: %s", lib_path))
    }
  if(!is.null(lib_file)) {
    if(!file.exists(file.path(lib_path, lib_file))) {
      stop(sprintf("\nFile not found: %s", file.path(lib_path, lib_file)))
    }
  }

  # read in blast result--------------------------------------------------------
  blast <- read.csv(blast_file, header = TRUE)
  blast <- blast[, c('query_num',	'sample_name',	'seq_len',
                     'hit_num',	'match',	'accession',	'total_score',
                     'perc_cover',	'perc_ident',	'evalue',
                     'match_description',	'query_seq')]

  # add time stamp
  blast <- cbind(date = Sys.Date(), blast)

  # if strain library not supplied, make new, empty library---------------------
  if(is.null(lib_file)) {

    lib_file <- file.path(lib_path, paste0(lib_name, ".csv"))

    # empty library df
    empty_lib <- blast[0,]

    write.csv(empty_lib, lib_file, row.names = FALSE, na = "")
    msg <- sprintf("No strain library supplied. Creating new strain library: \n%s\n",
                   lib_file)
    message(msg)
  }
  else {
    lib_file <- file.path(lib_path, lib_file)
  }

  # read in library------------------------------------------------------------
  lib_df <- read.csv(lib_file, header = TRUE, stringsAsFactors = FALSE)

  # if date columm doesn't exist, add it
  lib_df <- tryCatch({
    lib_df <- lib_df[, c('date','query_num',	'sample_name',	'seq_len',
                         'hit_num', 'match',	'accession',	'total_score',
                         'perc_cover', 'perc_ident',	'evalue',
                         'match_description',	'query_seq')]},
    error = function(e) {
      lib_df <- cbind(date = '', lib_df)
      return(lib_df)
    })

  # identify single species-----------------------------------------------------
  # put current library and current blast results together
  full_df <- suppressWarnings(gtools::smartbind(lib_df, blast))

  species_tally <- as.data.frame(table(full_df$match))
  colnames(species_tally) <- c('match','Freq')
  species_tally$match <- as.character(species_tally$match)

  single <- species_tally[species_tally$Freq == 1,]
  add_singles <- !single$match %in% lib_df$match

  if(any(add_singles)) {

    add_singles <- single$match[add_singles]

    msg <- sprintf("%s found as unique species that is not currently in your strain library. Adding to strain library.", add_singles)
    for(i in msg) message(i)

    # appending entry to library
    entry <- blast[blast$match %in% add_singles,]
    lib_df <- suppressWarnings(gtools::smartbind(lib_df, entry))
    write.csv(lib_df, lib_file, row.names = FALSE)
  }
  else {
    message("No unique species added to strain library.")
  }

  # Now go through rest of sequences and create fasta files for each species----


  # getting entries from strain library
  wip_df <- lib_df[lib_df$match %in% species_tally$match[species_tally$Freq > 1]
                   & lib_df$match %in% unique(blast$match), ]

  # getting entries from blast result
  wip_df <- suppressWarnings(gtools::smartbind(wip_df, blast[!blast$match %in% add_singles,]))

  msg <- sprintf("\nMultiple strains of %s were found when going through current blast file and strain library. Generating fasta file for this species for sequence alignment", unique(wip_df$match))

  for(i in msg) message(i)

  species_df <- wip_df[,c('sample_name','match','query_seq')] %>%
    group_by(match) %>%
    do(dummy = df2fasta(., fname = unique(.$match), path = lib_path,
                header_col = 'sample_name', seq_col = 'query_seq'))

}
