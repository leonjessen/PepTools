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

# Set a data tibble with names and abbreviations of the standard 20
# proteogenic amino acid residues
.set_aminoacids = function(){
  AMINOACIDS = matrix(c("Alanine","Ala","A","Arginine","Arg","R",
                        "Asparagine","Asn","N","Aspartate","Asp","D",
                        "Cysteine","Cys","C","Glutamine","Gln","Q",
                        "Glutamate","Glu","E","Glycine","Gly","G",
                        "Histidine","His","H","Isoleucine","Ile","I",
                        "Leucine","Leu","L","Lysine","Lys","K",
                        "Methionine","Met","M","Phenylalanine","Phe","F",
                        "Proline","Pro","P","Serine","Ser","S",
                        "Threonine","Thr","T","Tryptophan","Trp","W",
                        "Tyrosine","Tyr","Y","Valine","Val","V"),
                      ncol = 3, byrow = TRUE)
  colnames(AMINOACIDS) = c("full", "three", "one")
  AMINOACIDS = AMINOACIDS %>% as_tibble
  save(AMINOACIDS, file = "data/AMINOACIDS.RData")
  return(0)
}
.set_peptides = function(){
  set.seed(509279)
  p_file = paste0("https://raw.githubusercontent.com/leonjessen/",
                  "keras_tensorflow_demo/master/data/",
                  "ran_peps_netMHCpan40_predicted_A0201_reduced_cleaned_balanced.tsv")
  PEPTIDES = read_tsv(file = p_file)
  PEPTIDES = PEPTIDES %>%
    filter(label_chr=="SB") %>%
    sample_n(size = 5000) %>%
    pull(peptide)
  save(PEPTIDES, file = "data/PEPTIDES.RData")
  return(0)
}
# Download the blosum62 matrix from NCBI
# https://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
# and save it into data/BLOSUM62.RData
.get_blosum62 = function(){
  ncbi_url = "https://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt"
  BLOSUM62 = ncbi_url %>% read.table
  BLOSUM62 = BLOSUM62[1:20,1:20] %>% as.matrix
  save(BLOSUM62, file = "data/BLOSUM62.RData")
  return(0)
}
# Get the background frequencies
# Background distribution frequencies for Kullback-Leibler calculation is from:
#     Proteome-pI: proteome isoelectric point database
#     Lukasz P. Kozlowski
#     Nucleic Acids Research, Volume 45, Issue D1, 4 January 2017, Pages D1112–
#     D1116, https://doi.org/10.1093/nar/gkw978
#     Table 2. Amino acid frequency for the kingdoms of life in the Proteome-pI database
#     Ala 	Cys 	Asp 	Glu 	Phe 	Gly 	His 	Ile 	Lys 	Leu 	Met 	Asn
#     A     C     D     E     F     G     H     I     K     L     M     N
#     8.76 	1.38 	5.49 	6.32 	3.87 	7.03 	2.26 	5.49 	5.19 	9.68 	2.32 	3.93
#     Pro 	Gln 	Arg 	Ser 	Thr 	Val 	Trp 	Tyr   NA    NA
#     P     Q     R     S     T     V     W     Y     X     -
#     5.02 	3.90 	5.78 	7.14 	5.53 	6.73 	1.25 	2.91  4.55  4.55
#     'X' and '-' are set to 1/22*100 (i.e. flat bg)
# Table 2 large:
# https://academic.oup.com/view-large/51205195/Proteome-pI%3A%20proteome%20isoelectric%20point%20database
.get_bgfreqs = function(){
  BGFREQS = read_tsv(file = 'data/bgfreqs.txt', comment = '#')
  BGFREQS = BGFREQS %>% filter(Kingdom == 'All') %>%
    select(2:21) %>% as.matrix %>% .[1,]
  names(BGFREQS) = BGFREQS %>% names %>% residue_name
  BGFREQS = BGFREQS[AMINOACIDS$one] / 100
  save(BGFREQS, file = "data/BGFREQS.RData")
  return(0)
}
