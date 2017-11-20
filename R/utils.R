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



# ------------------------------------------------------------------------------
# COMPLETED FUNCTIONS
# ------------------------------------------------------------------------------
# Convert a vector of peptides to matrix form
pep_mat = function(pep){
  # Check input
  pep_check(pep = pep)
  # Convert to matrix
  # do.call applies a function to the list returned from args
  # so rbind to form matrix each of the elements in the list returned
  # by strsplit
  return( do.call(what = rbind, args = strsplit(x = pep, split = '')) )
}

# Generate 'n' random 'k'-mers (peptides)
pep_ran = function(n = 10, k = 9){

  # A bit tricky to avoid looping over millions of peptides

  # Define the standard 20 amino acid residue characters
  res_chars = AMINOACIDS$one

  # First we generate one long vector with all the samples residues
  smpl_chars = res_chars %>% sample(size = n*k, replace = TRUE)

  # Then we collapse into one long string
  smpl_string = smpl_chars %>% str_c(collapse = '')

  # Now we generate indices corresponding to extracting every 9 characters
  to_index   = seq(k, k*n, by = k)
  from_index = to_index - (k-1)

  # and extract using sub string
  kmer_peptides = smpl_string %>% str_sub(start = from_index, end = to_index)

  # Done! Return!
  return(kmer_peptides)
}

# Check if a vector of peptides is valid, according to criteria
#  - Is a vector with one or more elements
#  - Is a character vector
#  - All elements in vector have the same number of characters
#  - Elements only contain allowed amino acid residues 'ARNDCQEGHILKMFPSTWYVX-'
pep_check = function(pep){
  # Check if 'pep' is a character vector
  if( !( is.vector(pep) & is.character(pep) ) ){
    stop("'pep' has to be a character vector with one or more peptides!")
  }
  # Check if all sequences are of equal length
  if( pep %>% nchar %>% unique %>% length > 1 ){
    stop("All peptides must have equal number of positions")
  }
  # Check if 'pep' contain non-allowed amino acid residue characters
  if( pep %>% str_c(collapse='') %>% str_detect(pattern = "[^ARNDCQEGHILKMFPSTWYVX-]")){
    stop("Non standard amino acid residue character found.
         Only 'ARNDCQEGHILKMFPSTWYVX-' allowed")
  }
  return(TRUE)
}
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
# Score a peptide against a PSSM using sum
pep_score = function(pep, PSSM){
  pep_check(pep = pep)
  pssm_check(PSSM = PSSM)
  n_pep = length(pep)
  l_pep = nchar(pep)[1]
  if( l_pep != nrow(PSSM) ){
    stop("Number of positions in peptide must match number of rows in PSSM")
  }

  # The following is tricky in order to avoid looping over millions of peptides

  # Select columns from PSSM according to all peptides concatenated
  # to one long vector, where each element is a single residue
  pep_long  = pep %>% str_c(collapse = '') %>% str_split('') %>% unlist
  pssm_long = PSSM[,pep_long]

  # Now create repeated column bound diagonal matrices and multiply with
  # long PSSM effectively setting all non-diagonal scores to zero
  # Sum of columns will now be the score of each residue at each position
  s_vec = colSums( pssm_long * matrix(rep(diag(l_pep),n_pep), nrow = l_pep) )

  # Convert to matrix per l_pep positions
  s_mat = matrix(s_vec, ncol = l_pep, byrow = TRUE)

  # and compute pep sum scores as sum of rows
  scores = rowSums(s_mat)

  # Done
  return(scores)
}

# Split a peptide into residue characters
pep_split = function(pep){
  pep_check(pep = pep)
  if( length(pep) > 1 ){
    stop("Only one peptide is allowed. For mulitple peptides use pep_mat()")
  }
  return( pep %>% str_split('') %>% unlist )
}
# Mutate residues along peptide positions
pep_mutate = function(pep, pos, res){
  if( length(pep) > 1 ){
    stop("Only mutate one peptide at a time")
  }
  if( max(pos) > nchar(pep) ){
    stop("Trying to mutate positions beyond peptide length")
  }
  pep_check(pep = pep)
  p_mat = pep_mat(pep = pep)
  res = res %>% str_c(collapse='') %>% str_split('') %>% unlist
  mutants = sapply(1:nrow(p_mat),function(i){
    sapply(1:length(pos), function(j){
      sapply(1:length(res), function(k){
        mut_ijk = p_mat[i,]
        mut_ijk[pos[j]] = res[k]
        mut_ijk = mut_ijk %>% str_c(collapse = '')
        return(mut_ijk)
      })
    })
  })
  mutants = mutants %>% as.vector
  mutants = mutants[mutants != pep]
  return(mutants)
}
# Translate between amino acid residue symbold one <-> three <-> full
residue_name = function(res, to = 'one'){

  # Check if a tibble column was supplied
  if( is_tibble(res) ){
    res = res %>% pull
  }

  # Check arguments
  if( !any( to %in% c('one','three','full') )){
    stop("Argument 'to' must be one of:
         'one', 'three' or 'full'")
  }

  # Format input
  res  = res %>% tolower
  init = res %>% str_sub(start = 1, end = 1) %>% toupper
  res  = res %>% str_replace(pattern = "^\\w{1}", replacement = init)
  to   = to %>% tolower

  # Handle NAs
  nas_pos = c()
  if( any( is.na(res) ) ){
    nas_pos = which( is.na(res) )
    res[nas_pos] = sample(res[!is.na(res)],1) # Replace NAs with random element
  }

  # Hardcode tibble of amino acid residue names
  res_name_dat = AMINOACIDS

  # Check valid amino acid residue name in input
  if( !(all(res %in% res_name_dat$one) |
        all(res %in% res_name_dat$three) |
        all(res %in% res_name_dat$full)) ){
    stop("Non standard amino acid residue character or mix of found.
         Only one/three/full names for 'ARNDCQEGHILKMFPSTWYV' allowed")
  }

  # First we detect the input type supplied by the user
  input_type = res_name_dat %>% apply(2,function(x_j){ all(res %in% x_j) })
  input_type = colnames(res_name_dat)[input_type]

  # Then we define the names vector for translation
  trans_vec        = res_name_dat %>% pull(to)
  names(trans_vec) = res_name_dat %>% pull(input_type)

  # Return
  out = unname(trans_vec[res])
  out[nas_pos] = NA # Re-introduce replaced NAs
  return(out)

}
# Create empty matrix with corresponding row and column names
pssm_empty = function(npos = 9, res_chars = 'ARNDCQEGHILKMFPSTWYV', special_chars = FALSE){
  if( special_chars ){
    res_chars = paste0(res_chars, 'X-')
  }
  res_chars      = res_chars %>% str_split('') %>% unlist
  pssm           = matrix(data = 0, nrow = npos, ncol = length(res_chars))
  rownames(pssm) = seq(1, npos)
  colnames(pssm) = res_chars
  return(pssm)
}
# Scale a matrix, such that all values v lies in [0;1]
mat_zo_scale = function(x, na_rm = FALSE){
  if( !(is.matrix(x) & is.numeric(x)) ){
    stop("Input must be a numeric matrix")
  }
  x_scaled = ( x - min(x,na.rm=na_rm) ) / ( max(x,na.rm=na_rm) - min(x,na.rm=na_rm) )
  return(x_scaled)
}
# Scale a vector, such that all values v lies in [0;1]
vec_zo_scale = function(x, na_rm = FALSE){
  if( !(is.vector(x) & is.numeric(x)) ){
    stop("Input must be a numeric matrix")
  }
  x_scaled = ( x - min(x,na.rm=na_rm) ) / ( max(x,na.rm=na_rm) - min(x,na.rm=na_rm) )
  return(x_scaled)
}
# Translate all non-standard amino acid residue characters with 'X'
pep_clean = function(pep){
  return( pep %>% toupper %>%
            str_replace_all(pattern = "[^ARNDCQEGHILKMFPSTWYVX-]", replacement = "X") )
}
# Encode a vector of peptides using the BLOSUM62 matrix.
# Effectively each peptide become a 2D image and final
# encoding a 3D tensor
pep_encode = function(pep){

  # Check input vector
  pep_check(pep = pep)

  # First we retrive the BLOSUM62 matrix from data() and scale,
  # such that all values fall in the range 0 to 1
  blosum62 = BLOSUM62 %>% mat_zo_scale

  # Then we convert the vector of peptides to a matrix
  # with dimensions 'm x n' = n_peps x length_peps'
  p_mat = pep %>% pep_mat

  # Assign meaningful variable names to dimensions
  n_peps = length(pep)    # i'th row
  l_peps = nchar(pep[1])  # j'th column
  l_enc  = ncol(blosum62) # k'th slice

  # Finally we define our output tensor as a 3d array
  # with dimensions l_peps x l_enc x n_peps
  o_tensor = array(data = NA, dim = c(n_peps,l_peps, l_enc))
  for( i in 1:n_peps ){
    pep_i_residues = p_mat[i,]
    o_tensor[i,,]  = blosum62[pep_i_residues,]
  }
  return(o_tensor)
}
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
shannon_entropy = function(pep){
  pep_check(pep)
  f_mat = res_freqs(pep = pep)
  s_mat = -1 * f_mat * log2(f_mat)
  return(s_mat)
}

# As 4.2.1 but summed per position
shannon_entropy_pos = function(pep){
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

# From a list of peptides, calculate the corresponding frequency matrix
res_freqs = function(pep){
  pep_check(pep)
  c_mat = res_counts(pep = pep)
  f_mat = c_mat / rowSums(c_mat)
  return(f_mat)
}

# From a list of peptides, calculate the corresponding counts matrix
res_counts = function(pep){
  p_mat = pep_mat(pep)
  c_mat = pssm_empty(npos = ncol(p_mat))
  for( res in colnames(c_mat) ){
    c_mat[,res] = colSums(p_mat == res)
  }
  return(c_mat)
}
