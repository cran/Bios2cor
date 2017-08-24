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
rotamer_entropy <- function(rotamers){

  nb_rotamers <- length(rotamers[1,])
  nb_frames <- length(rotamers[,1])
  max_changes <- nb_frames-1
  res_names <- colnames(rotamers)

  entropy <- c()
  for(i in 1:nb_rotamers){
      changes <- 0
      
      for(j in 1:(nb_frames-1)){
	  if(rotamers[j,i] != rotamers[j+1, i]) changes <- changes+1
      }
      entropy[i] <- changes/max_changes
  }
  names(entropy) <- res_names
  return (entropy)
}

