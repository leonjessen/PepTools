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


