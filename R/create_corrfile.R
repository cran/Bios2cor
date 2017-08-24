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
create_corrfile <- function(corr_struct, entropy= NULL, auto_pairing= FALSE, filepath){
  
  ext <- strsplit(basename(filepath), "\\.")[[1]][2]
  filename <- strsplit(basename(filepath), "\\.")[[1]][1]
  filepath_auto <- paste(dirname(filepath), "/", filename, "_AUTO.", ext, sep= "")
  
  nb_pos <- length(corr_struct[1,])
  elem_names <- colnames(corr_struct)
  
  scale_struct <- (corr_struct-mean(corr_struct))/sd(corr_struct)
  
  angles_num <- c()
  if(auto_pairing == FALSE){
    ballesteros <- grepl("chi", elem_names, fixed= TRUE)
    ballesteros <- unlist(lapply(1:length(elem_names), function(i){if(ballesteros[i] == TRUE) 1 else 0}))
    nb_ballesteros <- sum(ballesteros)
    if(nb_ballesteros != 0){
      angles_num <- unlist(lapply(elem_names, function(x){strsplit(x, "[.]")[[1]][1]}))
    } else {
      angles_num <- elem_names
    }
  }
  
  head <- paste("pos_i", "pos_j", "entropy_i", "entropy_j", "corr", "z_score")
  write(head, file= filepath, append= FALSE)
  
  if(auto_pairing == TRUE){
    for(i in 1:nb_pos){
      for(j in i:nb_pos){
	pos_i <- elem_names[i]
	pos_j <- elem_names[j]
	if(pos_i != pos_j){
	  entropy_i <- entropy[pos_i]
	  entropy_j <- entropy[pos_j]
	  corr <- corr_struct[pos_i, pos_j]
	  z_score <- scale_struct[pos_i, pos_j]
	  
	  #Creating the line to insert to the file
	  current_line <- paste(pos_i, pos_j, entropy_i, entropy_j, corr, z_score)
	  
	  write(current_line, file= filepath_auto, append= TRUE)
	}
      }
    }
  } else {
    for(i in 1:nb_pos){
      for(j in i:nb_pos){
	ang_i <- angles_num[i]
	ang_j <- angles_num[j]
	if(ang_i != ang_j){
	  pos_i <- elem_names[i]
	  pos_j <- elem_names[j]
	  entropy_i <- entropy[pos_i]
	  entropy_j <- entropy[pos_j]
	  corr <- corr_struct[pos_i, pos_j]
	  z_score <- scale_struct[pos_i, pos_j]
	  
	  #Creating the line to insert to the file
	  current_line <- paste(pos_i, pos_j, entropy_i, entropy_j, corr, z_score)
	  
	  write(current_line, file= filepath, append= TRUE)
	}
      }
    }
  }
}

 
