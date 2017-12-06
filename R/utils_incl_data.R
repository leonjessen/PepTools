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
