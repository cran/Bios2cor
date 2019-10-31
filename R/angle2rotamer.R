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

angle2rotamer <- function(dynamic_structure, conversion_file=system.file("rotamer/dynameomics_rotamers.csv", package= "Bios2cor")){

if (missing(dynamic_structure)){
  stop("A 'dynamic_structure' object is required")
}


  #Importing rotamer information
  angles <- read.table(conversion_file, header= TRUE, fill= TRUE, sep= ",")

  #Importing torsional data from dynamic_structure object
  pdb <- dynamic_structure$pdb
  trj <- dynamic_structure$trj
  
  ca.inds <- dynamic_structure$ca.inds
  
  xyz <- dynamic_structure$xyz
  tor <- dynamic_structure$tor
  nb_torsions <- dynamic_structure$nb_torsions
  nb_frames <- dynamic_structure$nb_frames
  
  prot.seq <- dynamic_structure$prot.seq
  
  tor.names <- dynamic_structure$tor.names
  tor.resno <- dynamic_structure$tor.resno
  tor.angle <- dynamic_structure$tor.angle
  tor.seq <- dynamic_structure$tor.seq

  converted_angles <- matrix("", ncol= nb_torsions, nrow= nb_frames)
  colnames(converted_angles) <- tor.names
  
  #converting dihedral angles into rotamers 
  for(i in 1:nb_frames){
       for(j in 1:nb_torsions){
	  residue <- tor.seq[j]
	  dihedral <- tor.angle[j]
	  angle <- tor[i,j]
	  
	  #select the rotamers possible for the dihedral angle and residue under consideration
	  selected_torsion <- as.matrix(angles[(angles[,1] == residue) & (angles[,2] == dihedral),])
	  rotamer_nb <- length(selected_torsion[,1])
	  
	  #check if the selected torsion is contained in the conversion file
	  if(rotamer_nb != 0){
	    for(r in 1:rotamer_nb){
		ang_1 <- selected_torsion[r,4]
		ang_2 <- selected_torsion[r,5]
		if(ang_1 == "0-" | ang_1 == "0+") ang_1 <- "0"
		if(ang_2 == "0-" | ang_2 == "0+") ang_2 <- "0"
		
		min_angle <- as.numeric(ang_1)
		max_angle <- as.numeric(ang_2)
		
		special_interval <- FALSE #the interval is greater than 180 degrees
		if(min_angle > max_angle) special_interval <- TRUE
		
		if(special_interval){
		    if((min_angle <= angle & angle <= 180) | (-180 <= angle & angle <= max_angle)) converted_angles[i,j] <- selected_torsion[r,3]
		} else {
		    if(min_angle <= angle & angle <= max_angle) converted_angles[i,j] <- selected_torsion[r,3]
		}
	    }
	  } else {
	      cat(paste("dihedral angle not found in the conversion file \n"))
	      cat(paste("residue = ", residue, "\n"))
	      cat(paste("dihedral = ", dihedral, "\n"))
	  }
      }
  }
  
  return (converted_angles)
}
