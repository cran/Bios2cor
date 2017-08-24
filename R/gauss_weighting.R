# Package: Bios2cor 
# This file is part of Bios2cor R package.
# Bios2cor is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# See the GNU General Public License at:
# http://www.gnu.org/licenses/
#

#L is recommanded to be between 0.05 and 0.2, it corresponds to gaussian curve's width
#lambda must be between 0 and 1
gauss_weighting <- function(entropy, L= 0.1, lambda= 0){
  nb_pos <- length(entropy)
  res_names <- names(entropy)
   
  res <- matrix(1/nb_pos, nrow= nb_pos, ncol= 1)
  
  weighting <- unlist(lapply(1:nb_pos, function(pos){ exp(-(entropy[pos]-0.5)^2/L) }))
  
  #adjusting weighting
  weighting <- lambda + (1-lambda) * weighting

  #Normalization
  SUM <- sum(weighting)
  res <- weighting/SUM
  
  names(res) <- res_names
  return (res)
}
 
