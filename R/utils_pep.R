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
#' Convert a vector of peptides to matrix form
#'
#' @param pep A vector of equal length peptides
#' @return A n_pep x l_pep matrix
#' @examples
#' pep_mat(c("RQGQDHPTM","RGQKTTDNA","NILYEYWDY"))
pep_mat = function(pep){
  # Check input
  pep_check(pep = pep)
  # Convert to matrix
  # do.call applies a function to the list returned from args
  # so rbind to form matrix each of the elements in the list returned
  # by strsplit
  return( do.call(what = rbind, args = strsplit(x = pep, split = '')) )
}
################################################################################
################################################################################
################################################################################





################################################################################
################################################################################
################################################################################
#' Generate 'n' random 'k'-mers (peptides)
#'
#' @param n The number of peptide(s) to generate
#' @param k The length of the peptide(s)
#' @return A character vector of n peptides
#' @examples
#' pep_ran(n = 1000, k = 9)
pep_ran = function(n, k){

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
################################################################################
################################################################################
################################################################################





################################################################################
################################################################################
################################################################################
#' Check if a vector of peptides is valid.
#'
#' According to criteria:
#'  i. Is a vector with one or more elements.
#'  ii. Is a character vector.
#'  iii. All elements in vector have the same number of characters.
#'  iv. Elements only contain allowed amino acid residues \code{'ARNDCQEGHILKMFPSTWYVX-'}.
#'
#' @param pep A vector of peptides to be checked
#' @return TRUE if the vector passed the check, otherwise stop() is called
#' @examples
#' pep_check(c("RQGQDHPTM","RGQKTTDNA","NILYEYWDY"))
#' pep_check(c("RQGQDHPTM","RGQKTTDNA","NILYEYWDYX"))
#' pep_check(c("RQGQDHPTM","RGQKTTDNA","NILYEYWDB"))
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
################################################################################
################################################################################
################################################################################





################################################################################
################################################################################
################################################################################
#' Score a peptide against a PSSM using sum of positional scores.
#'
#' @param pep A vector of equal length peptides to be scored
#' @param PSSM A Position Specific Scoring Matrix, where the number of rows must be equal to the length of the input peptides.
#' @return A vector of scores, one for each input peptides. Scores are sorted according to input peptides.
#' @examples
#' PSSM = pssm_empty(npos = 9)
#' PSSM[1:9,1:20] = rnorm(180)
#' pep_score(c("RQGQDHPTM","RGQKTTDNA","NILYEYWDY"),PSSM)
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
################################################################################
################################################################################
################################################################################





################################################################################
################################################################################
################################################################################
#' Split a peptide into residue characters
#' @param pep A string representing a peptide
#' @return A character vector with number of elements corresponding to the number of amino acids in the input peptide
#' @examples
#' pep_split("RQGQDHPTM")
pep_split = function(pep){
  pep_check(pep = pep)
  if( length(pep) > 1 ){
    stop("Only one peptide is allowed. For mulitple peptides use pep_mat()")
  }
  return( pep %>% str_split('') %>% unlist )
}
################################################################################
################################################################################
################################################################################





################################################################################
################################################################################
################################################################################
#' Mutate residues along peptide positions
#' @param pep A string representing a peptide
#' @param pos Which position(s) to mutate
#' @param res Which amino acid(s) to mutate to
#' @param keep_wildtype Boolean should \code{pep} be included in the output
#' @return A character vector with number of elements corresponding to the
#' number of peptides required to cover all mutations
#' @examples
#' pep_mutate(pep = "RQGQDHPTM", pos = seq(1,9), res = AMINOACIDS$one)
pep_mutate = function(pep, pos, res, keep_wt = FALSE){
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
  mutants = mutants %>% as.vector %>% unique
  if( !keep_wt ){
    mutants = mutants[mutants != pep]
  }
  if( length(mutants) == 0 ){
    stop("No peptides could be returned, e.g. pep_mutate('AAA',1:3,'A',keep_wt=FALSE)")
  }
  return(mutants)
}
################################################################################
################################################################################
################################################################################





################################################################################
################################################################################
################################################################################
#' Translate between amino acid residue symbols
#' @param res A character vector of residue symbols either in one-, three- or
#' full format
#' @param to Which symbol to translate to, one of \code{'one', 'three', 'full'}
#' @return A character vector with the translated symbols
#' @examples
#' residue_name(AMINOACIDS$full, to = 'three')
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
################################################################################
################################################################################
################################################################################





################################################################################
################################################################################
################################################################################
#' Translate all non-standard amino acid residue characters with 'X'
#' @param pep A character vector of peptides to be cleaned
#' @return A character vector of cleaned peptides
#' @examples
#' pep_clean(sample(letters,20,replace=TRUE))
pep_clean = function(pep){
  return( pep %>% toupper %>%
            str_replace_all(pattern = "[^ARNDCQEGHILKMFPSTWYVX-]", replacement = "X") )
}
################################################################################
################################################################################
################################################################################





################################################################################
################################################################################
################################################################################
#' Encode a peptide
#'
#' Peptides are encoded using the BLOSUM62 probability matrix
#' (See BLOSUM62_PROB data matrix)
#'
#' Each position in the peptide become a vector of 20 values, corresponding to
#' the rounded log odds ratio for substituting the amino acid in the peptide
#' with each of the 20 standard proteogenic amino acids.
#'
#' The final result is a 3D array (tensor) of peptide 'images' with \code{n} rows, \code{k}
#' columns and \code{20} slices, which is stored as a 3D array (tensor)
#'
#' @param pep A character vector of peptides to be encoded
#' @return A 3D array (tensor) of peptide 'images'
#' @examples
#' pep_encode(pep_ran(k=9,n=10))
#' dim(pep_encode(pep_ran(k=9,n=10)))
pep_encode = function(pep){

  # Check input vector
  pep_check(pep = pep)

  # Set encoding matrix
  bl62_prob = BLOSUM62_PROB

  # Then we convert the vector of peptides to a matrix
  # with dimensions 'm x n' = 'n_peps x length_peps'
  p_mat = pep %>% pep_mat

  # Assign meaningful variable names to dimensions
  n_peps = length(pep)    # i'th row
  l_peps = nchar(pep[1])  # j'th column
  l_enc  = ncol(bl62_prob) # k'th slice

  # Finally we define our output tensor as a 3d array
  # with dimensions n_peps x l_peps x l_enc (l_enc = 20)
  o_tensor = array(data = NA, dim = c(n_peps,l_peps, l_enc))
  for( i in 1:n_peps ){
    pep_i_residues = p_mat[i,]
    pep_img = bl62_prob[pep_i_residues,]
    o_tensor[i,,]  = pep_img
  }
  return(o_tensor)
}
################################################################################
################################################################################
################################################################################





################################################################################
################################################################################
################################################################################
#' Get peptide images
#'
#' Peptides are encoded using \code{pep_encode()}, resulting in peptide 'images'
#'
#' Each position in the peptide become a vector of 20 values, corresponding to
#' the rounded log odds ratio for substituting the amino acid in the peptide
#' with each of the 20 standard proteogenic amino acids. The final result is a
#' list of peptide 'images'
#'
#' @param pep A character vector of peptides to be converted to 'images'
#' @return A list of peptide 'images'
#' @examples
#' pep_get_images(pep_ran(k=9,n=10))
pep_get_images = function(pep){

  # Check input vector
  pep_check(pep = pep)

  # Set encoding matrix
  bl62_prob = BLOSUM62_PROB

  # Then we convert the vector of peptides to a matrix
  # with dimensions 'm x n' = 'n_peps x length_peps'
  p_mat = pep %>% pep_mat

  # Generate peptide 'images' and save in list
  o_list = list()
  for( i in 1:nrow(p_mat) ){
    pep_i_residues = p_mat[i,]
    pep_img = bl62_prob[pep_i_residues,]
    o_list[[i]] = pep_img
  }
  return(o_list)
}
################################################################################
################################################################################
################################################################################





################################################################################
################################################################################
################################################################################
#' Plot peptide images
#'
#' Plot the first \code{n} encoded peptide 'images'
#'
#' Peptides are encoded using \code{pep_encode()}, resulting in peptide 'images'
#' Each position in the peptide become a vector of 20 values, corresponding to
#' the rounded log odds ratio for substituting the amino acid in the peptide
#' with each of the 20 standard proteogenic amino acids. The final result is a
#' list of peptide 'images' each with values in \code{[0;1]}.
#'
#' @param pep A character vector of peptides to be plotted as 'images'
#' @return A cowplot plot grid of the \code{n} first peptide 'images'
#' @examples
#' pep_plot_images(pep_ran(k=9,n=10))
pep_plot_images = function(pep, n = 3){

  # Check input vector
  pep_check(pep = pep)

  # Check 'n'
  if( n > length(pep) ){
    n = length(pep)
  }
  if( n < 1 ){
    stop("'n' must be larger than or equal to 1")
  }

  # Convert peptide to list of 'images'
  pep_imgs = pep[1:n] %>% pep_get_images

  # Convert each image in list to a plot
  plot_list = list()
  for( i in 1:length(pep_imgs) ){
    pep_img  = pep_imgs[[i]]
    residues = paste(rownames(pep_img), seq(1,nrow(pep_img)), sep = '_')
    residues = factor(x = residues, levels = rev(residues))
    p_img = tibble(pep_res = residues) %>%
      bind_cols(as_tibble(pep_img)) %>%
      gather(sub_res,val,-pep_res) %>%
      select(pep_res,sub_res,val) %>%
      ggplot(aes(sub_res,pep_res)) +
      geom_tile(aes(fill = val), color = "white") +
      scale_fill_gradient(low = "white", high = "black", limits = c(0,1)) +
      xlab("Score for substituting peptide residue with ...") +
      ylab("Peptide residue") +
      theme_classic() +
      theme(axis.line.x = element_blank(),
            axis.line.y = element_blank()) +
      ggtitle(paste("Peptide",pep[i],"'image' encoded"))
    plot_list[[i]] = p_img
  }

  # Done return cowplot plotlist
  return( plot_grid(plotlist = plot_list, ncol = 1) )

}
################################################################################
################################################################################
################################################################################





################################################################################
################################################################################
################################################################################
#' Retrieve the consensus sequence from a set of peptides
#' @param pep A vector of peptides
#' @return A string corresponding to the consensus sequence
#' @examples
#' pep_consensus(c("YMNSMQEML","FIYRHMFCV","VLFKFDMFI","KLLDRFPVA","RVLDDFTKL"))
pep_consensus = function(pep){
  get_max = function(x){
    o = names(x)[which.max(x)]
    if( length(o) > 1 ){
      o = paste(o,collapse='')
      o = paste0("[",o,"]")
    }
    return(names(x)[which.max(x)])
  }
  pep_check(pep = pep)
  pep %>% pep_mat %>% apply(2,function(p_i){ p_i %>% table %>% return }) %>%
    lapply(get_max) %>% unlist %>% paste(collapse = '') %>% return
}
################################################################################
################################################################################
################################################################################





################################################################################
################################################################################
################################################################################
#' Encode a peptide
#'
#' Peptides are encoded using the BLOSUM62 probability matrix
#' (See BLOSUM62_PROB data matrix)
#'
#' Each position in the peptide become a vector of 20 values, corresponding to
#' the rounded log odds ratio for substituting the amino acid in the peptide
#' with each of the 20 standard proteogenic amino acids.
#'
#' The final result is a tibble with as many rows as input peptides and number
#' columns equal to the number of positions in the input peptides times the
#' number of BLOSUM62 substitions values (20) plus a column with the input
#' peptides. For 100 9-mers this yields 100 rows and 181 columns
#'
#' @param pep A character vector of peptides to be encoded
#' @return A tibble of encoded peptides
#' @examples
#' pep_encode_mat(pep_ran(k=9,n=100))
#' dim(pep_encode_mat(pep_ran(k=9,n=100)))
pep_encode_mat = function(pep){

  # Check input
  pep_check(pep)

  # Get number of positions in peptide
  n_pos = nchar(pep[1])

  # Set encoding matrix
  bl62_prob = BLOSUM62_PROB

  # Convert input peptides to one long character vector of single amino acids
  pep_str = pep %>% paste(collapse = '') %>% str_split('') %>% unlist

  # Generate output matrix, such that each position in the peptide is encoded as
  # a vector of 20 values corresponding to the values in the scaled BLOSUM62
  # matrix. Each position is then concatenated forming a vector of a length of
  # of number of positions in the peptides times the number of substitution
  # values in the BLOSUM62 matrix, i.e. for a 9-mer this will be 9 x 20 = 180.
  # The output matrix will then have dimensions 180 columns times number of input
  # peptides rows and then an extra column when includeing the input peptides
  # in the output. This way column names are:
  # A_p1 R_p1 N_p1 D_p1 C_p1 Q_p1 E_p1 G_p1 H_p1 I_p1 L_p1 K_p1 M_p1 F_p1 P_p1
  # S_p1 T_p1  W_p1  Y_p1  V_p1  A_p2 R_p2 N_p2 D_p2 C_p2 Q_p2 E_p2 G_p2 ...
  out_mat = bl62_prob[pep_str,] %>% t %>% as.vector %>%
    matrix(ncol=n_pos * ncol(bl62_prob), byrow = TRUE)
  colnames(out_mat) = paste(colnames(bl62_prob),
                            paste0('p',rep(1:n_pos,rep(ncol(bl62_prob),n_pos))),sep='_')
  out_mat = out_mat %>% as_tibble
  out_mat = as_tibble(pep) %>% bind_cols(out_mat) %>% rename(peptide=value)
  return(out_mat)
}
################################################################################
################################################################################
################################################################################
