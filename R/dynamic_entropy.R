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
dynamic_entropy <- function(rotamers){

  if(missing(rotamers)){
     stop("A 'rotamers' matrix is required")
  }
 
  if(!is.matrix(rotamers)){
     stop("The argument is not a matrix")
  }

  nb_rotamers <- length(rotamers[1,])
  nb_frames <- length(rotamers[,1])
  max_changes <- nb_frames-1
  dihed_names <- colnames(rotamers)

  entropy <- matrix(0, ncol= nb_rotamers, nrow= 1)

  for(i in 1:nb_rotamers){
      changes <- 0
      for(j in 1:(nb_frames-1)){
	  if(rotamers[j,i] != rotamers[j+1, i]) changes <- changes+1
      }
      entropy[i] <- changes/max_changes
  }
  names(entropy) <- dihed_names
  return (entropy)
}

