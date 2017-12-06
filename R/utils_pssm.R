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
# Check if a PSSM has the correct format, according to criteria:
#  - Is a matrix
#  - Is numeric
#  - Has 20, 21 or 22 columns
#  - Column names only contain allowed amino acid residues 'ARNDCQEGHILKMFPSTWYVX-'
pssm_check = function(PSSM){
  if( !is.matrix(PSSM) | !is.numeric(PSSM) ){
    stop("PSSM must be a numeric matrix")
  }
  if( ncol(PSSM) < 20 | ncol(PSSM) > 22 ){
    stop("Number of columns in PSSM should be 20, 21 or 22:
         'ARNDCQEGHILKMFPSTWYV', 'ARNDCQEGHILKMFPSTWYV-' or 'ARNDCQEGHILKMFPSTWYVX-'
         Format: Rows = Positions in peptide, columns = residue characters")
  }
  if( PSSM %>% colnames %>% str_c(collapse='') %>% str_detect("[^ARNDCQEGHILKMFPSTWYVX-]") ){
    stop("Invalid column names in PSSM. Only 'ARNDCQEGHILKMFPSTWYVX-' allowed")
  }
  return(TRUE)
  }
# From a matrix of frequencies, where the peptide positions are the rows and
# where the amino acid residues are the columns and the sum of the rows all
# are equal to one, calculate the weighted bits of information
pssm_bits = function(PSSM){

  # Check that the input PSSM is a numeric matrix
  if( !is.matrix(PSSM) | !is.numeric(PSSM) ){
    stop("PSSM must be a numeric matrix")
  }

  # Check that the input PSSM is a frequency matrix
  if( !all(PSSM %>% rowSums() %>% round(2) == 1) ){
    stop("PSSM must be a frequency matrix, where all row sums are 1")
  }

  # Compute the postional relative Entropy (information content)
  ic = log2(20) + rowSums( PSSM * log2(PSSM) , na.rm = TRUE)

  # Compute bits of information for all PSSM[i,j]
  bits = PSSM * ic

  # check
  if( !all( round(ic,2) == round(rowSums(bits),2) ) ){
    stop("There was an error in computing the bits!")
  }

  # Done, return!
  return(bits)
}
# Create empty matrix with corresponding row and column names
pssm_empty = function(npos = 9){
  res_chars      = AMINOACIDS$one
  pssm           = matrix(data = 0, nrow = npos, ncol = length(res_chars))
  rownames(pssm) = seq(1, npos)
  colnames(pssm) = res_chars
  return(pssm)
}
# Create empty matrix with corresponding row and column names
pssm_empty_long = function(npos = 9){
  res_chars      = c(AMINOACIDS$one,'X','-')
  pssm           = matrix(data = 0, nrow = npos, ncol = length(res_chars))
  rownames(pssm) = seq(1, npos)
  colnames(pssm) = res_chars
  return(pssm)
}
# Flatten PSSM to a single row, such that 9 x 20 => 1 x 180
pssm_flatten = function(PSSM){
  out_col_names = sapply(1:ncol(PSSM),function(i){
    paste(rownames(PSSM),colnames(PSSM)[i],sep="_")
  }) %>% t %>% matrix(nrow=1)
  out = PSSM %>% t %>% matrix(nrow = 1)
  colnames(out) = out_col_names
  return(out)
}
# From a list of peptides, calculate the corresponding frequency matrix
pssm_freqs = function(pep){
  pep_check(pep)
  c_mat = pssm_counts(pep = pep)
  f_mat = c_mat / rowSums(c_mat)
  return(f_mat)
}
# From a list of peptides, calculate the corresponding counts matrix
pssm_counts = function(pep){
  p_mat = pep_mat(pep)
  c_mat = pssm_empty(npos = ncol(p_mat))
  for( res in colnames(c_mat) ){
    c_mat[,res] = colSums(p_mat == res)
  }
  return(c_mat)
}
# ------------------------------------------------------------------------------
# WORK IN PROGRESS
# ------------------------------------------------------------------------------
.correlate_pssms = function(PSSM_list, n = 10000, k = 9){
  if( !is_list(PSSM_list) ){
    stop("Input to correlate_pssms must be a list")
  }
  test_peps = pep_ran(n = n, k = k)
  out = matrix(test_peps,ncol=1)
  ids = names(PSSM_list)
  for( id in ids ){
    pssm = PSSM_list[[id]]
    out = cbind(out, pep_score(pep = test_peps, PSSM = pssm))
  }
  out = out[,-1]
  colnames(out) = ids
  out = apply(out,2,function(x_j){ return(as.numeric(x_j)) })
  rownames(out) = test_peps

  return(out)
}

# Information metric equations are taken from
#     O. Lund, S. Brunak, C. Kesmir, C. Lundegaard, M. Nielsen. Immunological
#     Bioinformatics. The MIT Press, Cambridge, Massachusetts, London, England,
#     1 st edition, 2005. ISBN: 9780262122801
#     Chapter 4: Methods Applied in Immunological Bioinformatics
#     Section 4.2: Information Carried by Immunogenic Sequences
#
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

# 4.2.1 Entropy
# Eq. 4.7, p70
.shannon_entropy = function(pep){
  pep_check(pep)
  f_mat = res_freqs(pep = pep)
  s_mat = -1 * f_mat * log2(f_mat)
  return(s_mat)
}

# As 4.2.1 but summed per position
.shannon_entropy_pos = function(pep){
  pep_check(pep)
  s_mat = shannon_entropy(pep)
  s_vec = rowSums(s_mat)
  names(s_vec) = rownames(s_mat)
  return(s_vec)
}

# 4.2.2 Relative Entropy
# Eq. 4.8, p70
.kullback_leibler_divergence = function(x){
  bgfreq  = c('A'=8.76,'C'=1.38,'D'=5.49,'E'=6.32,'F'=3.87,
              'G'=7.03,'H'=2.26,'I'=5.49,'K'=5.19,'L'=9.68,
              'M'=2.32,'N'=3.93,'P'=5.02,'Q'=3.90,'R'=5.78,
              'S'=7.14,'T'=5.53,'V'=6.73,'W'=1.25,'Y'=2.91,
              'X'=4.55,'-'=4.55) / 100
  pep_mat = pep_mat(x = x)
  out_mat = pssm_empty(npos = ncol(pep_mat))
  for( i in 1:nrow(out_mat) ){
    p_i = pep_mat[,i]
    for( j in 1:ncol(out_mat) ){
      res  = colnames(out_mat)[j]
      prob = mean(p_i==res)
      if( prob > 0 ){
        q = bgfreq[res]
        out_mat[i,j] = prob * log2(prob/q) # D(p||q)
      }
    }
  }
  return(out_mat)
}

# As 4.2.2 but summed per position
.kullback_leibler_divergence_pos = function(x){
  kld_mat = kullback_leibler_divergence(x)
  kld_vec = apply(kld_mat,1,function(p_i){ sum(p_i) })
  names(kld_vec) = paste0('p',1:length(kld_vec))
  return(kld_vec)
}

# 4.2.3 Logo Visualization of Relative Entropy
# Eq. 4.8, p71
.information_content = function(x){
  pep_mat = pep_mat(x = x)
  out_mat = pssm_empty(npos = ncol(pep_mat))
  for( i in 1:nrow(out_mat) ){
    p_i = pep_mat[,i]
    for( j in 1:ncol(out_mat) ){
      res  = colnames(out_mat)[j]
      prob = mean(p_i==res)
      if( prob > 0 ){
        out_mat[i,j] = prob * log2(prob/(1/20)) # I
      }
    }
  }
  return(out_mat)
}

# As 4.2.3 but summed per position
.information_content_pos = function(x){
  ic_mat = information_content(x)
  ic_vec = apply(ic_mat,1,function(p_i){ sum(p_i) })
  names(ic_vec) = paste0('p',1:length(ic_vec))
  return(ic_vec)
}


