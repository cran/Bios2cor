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

angle_conversion <- function(dynamic_struct, conversion_file){

  angles <- read.table(conversion_file, header= TRUE, fill= TRUE, sep= ",")
  print(head(angles))

  #Importing rotamer informations
  pdb <- dynamic_struct$pdb
  trj <- dynamic_struct$trj
  
  ca.inds <- dynamic_struct$ca.inds
  
  xyz <- dynamic_struct$xyz
  nb_frames <- dynamic_struct$nb_frames
  
  tor <- dynamic_struct$tor
  nb_residus <- dynamic_struct$nb_residus
  nb_frames <- dynamic_struct$nb_frames
  
  prot.seq <- dynamic_struct$prot.seq
  
  tor.names <- dynamic_struct$tor.names
  tor.resno <- dynamic_struct$tor.resno
  tor.angle <- dynamic_struct$tor.angle
  tor.seq <- dynamic_struct$tor.seq

  converted_angles <- matrix("", ncol= nb_residus, nrow= nb_frames)
  colnames(converted_angles) <- tor.names
  
  for(i in 1:nb_frames){
      cat(paste("frame", i, "\n"))
      for(j in 1:nb_residus){
	  residu <- tor.seq[j]
	  chain <- tor.angle[j]
	  angle <- tor[i,j]
	  
	  #reducing the total matrix to frames that contain the current residue (e.g. "300") and chain (e.g. "chi1")
	  reduced_angles <- as.matrix(angles[(angles[,1] == residu) & (angles[,2] == chain),])
	  nb_red <- length(reduced_angles[,1])
	  
	  #check if the current residu is contained in the conversion file
	  if(nb_red != 0){
	    for(r in 1:nb_red){
		ang_1 <- reduced_angles[r,4]
		ang_2 <- reduced_angles[r,5]
		if(ang_1 == "0-" | ang_1 == "0+") ang_1 <- "0"
		if(ang_2 == "0-" | ang_2 == "0+") ang_2 <- "0"
		
		min_angle <- as.numeric(ang_1)
		max_angle <- as.numeric(ang_2)
		
		special_interval <- FALSE #the interval is greater than 180 degrees
		if(min_angle > max_angle) special_interval <- TRUE
		
		if(special_interval){
		    if((min_angle <= angle & angle <= 180) | (-180 <= angle & angle <= max_angle)) converted_angles[i,j] <- reduced_angles[r,3]
		} else {
		    if(min_angle <= angle & angle <= max_angle) converted_angles[i,j] <- reduced_angles[r,3]
		}
	    }
	  } else {
	      cat(paste("Residu not found\n"))
	      cat(paste("residu = ", residu, "\n"))
	      cat(paste("chain = ", chain, "\n"))
	  }
      }
  }
  
  return (converted_angles)
}
