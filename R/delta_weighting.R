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

delta_weighting <- function(entropy, Smin= 0.4, Smax= 0.6){
  nb_pos <- length(entropy)
  res_names <- names(entropy)
  
  res <- matrix(1/nb_pos, nrow= nb_pos, ncol= 1)
  
  weighting <- unlist(lapply(1:nb_pos, function(pos){ if(entropy[pos] <= Smax & entropy[pos] >= Smin) 1 else 0 }))
  
  #Used to avoid 0 values
  eps <- 0.1
  weighting <- weighting + eps
  
  #Normalization
  SUM <- sum(weighting)
  res <- weighting/SUM
  
  names(res) <- res_names
  return (res)
}
 
