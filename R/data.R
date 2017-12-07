#' Amino Acid Residue Background Frequencies for Kongdom All
#'
#' Background distribution frequencies for Kullback-Leibler calculation is from:
#'     Proteome-pI: proteome isoelectric point database
#'     Lukasz P. Kozlowski
#'     Nucleic Acids Research, Volume 45, Issue D1, 4 January 2017, Pages D1112–
#'     D1116, https://doi.org/10.1093/nar/gkw978
#'
#' @format A named numeric vector with 20 elements:
#' \describe{
#'   \item{name}{One-letter symbol for each of the 20 standard proteogenic amino acids}
#'   \item{value}{The background frequency for kingdom all}
#' }
#' @source \url{https://academic.oup.com/nar/article/45/D1/D1112/2333931}
"BGFREQS"

#' Symbol translation table for the 20 standard proteogenic amino acids
#'
#' Symbol translation table for the 20 standard proteogenic amino acids
#'
#' @format A tibble with 20 rows and 3 columns:
#' \describe{
#'   \item{full}{Full-name symbol}
#'   \item{three}{Three-letter symbol}
#'   \item{one}{One-letter symbol}
#' }
"AMINOACIDS"

#' 5,000 9-mer peptides
#'
#' 5,000 9-mer peptides predicted to be strong binders to HLA-A*02:01 by netMHCpan-4.0
#'
#' @format A character vector with 5,000 strings of length 9
"PEPTIDES"

#' The BLOSUM62 matrix
#'
#' The BLock SUbstition matrix for sequences with less than 62\% similarity
#'
#' @format A named numeric matrix with 20 rows and 20 columns
#' @source \url{https://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt}
"BLOSUM62"
