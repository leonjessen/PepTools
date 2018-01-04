################################################################################
# PepTools - An R-package for making immunoinformatics accessible              #
#     Copyright (C) 2017 Leon Eyrich Jessen                                    #
################################################################################
# This program is free software: you can redistribute it and/or modify         #
# it under the terms of the GNU General Public License as published by         #
# the Free Software Foundation, either version 3 of the License, or            #
# (at your option) any later version.                                          #
#                                                                              #
# This program is distributed in the hope that it will be useful,              #
# but WITHOUT ANY WARRANTY; without even the implied warranty of               #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                #
# GNU General Public License for more details.                                 #
#                                                                              #
# You should have received a copy of the GNU General Public License            #
# along with this program.  If not, see <https://www.gnu.org/licenses/>.       #
################################################################################





################################################################################
################################################################################
################################################################################
#' Min-max scale matrix
#'
#' Min-max scale matrix, such that all values v lie in [0;1]. If X is a numeric
#' matrix, then min-max scale is: ( X - min(X) ) / ( max(X) - min(X) )
#' @param x A numeric matrix
#' @param na_rm Boolean: Should NAs be removed
#' @examples
#' mat_mima_scale(matrix(rnorm(180),nrow=9,ncol=20))
mat_mima_scale = function(x, na_rm = FALSE){
  if( !(is.matrix(x) & is.numeric(x)) ){
    stop("Input must be a numeric matrix")
  }
  x_scaled = ( x - min(x,na.rm=na_rm) ) / ( max(x,na.rm=na_rm) - min(x,na.rm=na_rm) )
  return(x_scaled)
}
################################################################################
################################################################################
################################################################################





################################################################################
################################################################################
################################################################################
#' Min-max scale vector
#'
#' Min-max scale vector, such that all values v lie in [0;1]. If x is a numeric
#' vector, then min-max scale is: ( x - min(x) ) / ( max(x) - min(x) )
#' @param x A numeric vector
#' @param na_rm Boolean: Should NAs be removed
#' @examples
#' vec_mima_scale(rnorm(20))
vec_mima_scale = function(x, na_rm = FALSE){
  if( !(is.vector(x) & is.numeric(x)) ){
    stop("Input must be a numeric matrix")
  }
  x_scaled = ( x - min(x,na.rm=na_rm) ) / ( max(x,na.rm=na_rm) - min(x,na.rm=na_rm) )
  return(x_scaled)
}
################################################################################
################################################################################
################################################################################





################################################################################
################################################################################
################################################################################
#' Read a FASTA file
#'
#' read_fasta reads a FASTA file of nucleotides or amino acids file and returns
#' a tibble with number of rows corresponding to the number of sequences and two
#' variables: 'fasta_header' and 'sequences'
#' @param file The FASTA file to be read either local or URL
#' @examples
#' read_fasta(file = 'my_fasta_file.fsa')
#' read_fasta(file = 'https://www.ncbi.nlm.nih.gov/WebSub/html/help/sample_files/nucleotide-sample.txt')
read_fasta = function(file){

  # Check input
  if( !is_character(file) ){
    stop("'file' has to be a string specifying a file name")
  }
  #if( !file.exists(file) ){
  #  stop(paste("Unable to read file",file))
  #}

  # Read lines from file
  lines = file %>% read_lines

  # Remove any comments
  keep  = lines %>% str_detect('#') %>% !.
  lines = lines[keep]
  rm(keep)

  # Set variables
  n_seqs    = lines %>% str_detect("^>") %>% sum
  headers   = rep(NA, n_seqs)
  sequences = rep('', n_seqs)

  # Iterate over lines and fill containers
  entry_no = 1
  for( line in lines ){
    if( line %>% str_detect("^>") ){
      headers[entry_no] = line
      entry_no = entry_no + 1
    } else {
      sequences[entry_no-1] = paste0(sequences[entry_no-1], line)
    }
  }

  # Delete lines
  rm(lines)

  # Return output tibble
  return(tibble(fasta_header = headers, sequence = sequences))
}
################################################################################
################################################################################
################################################################################
