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

rotamer_omes <- function(dynamic_struct, rotamers, res_selection= c("C","I","L","M","V","R","H","K","D","E","N","Q","F","Y","W","T","S","P"), z_score= TRUE, auto_pairing= FALSE){
  nb_angles <- length(rotamers[1,])
  nb_frames <- length(rotamers[,1])
  
  pos_names <- colnames(rotamers)
  
  names <- unique(as.vector(rotamers))
  nb_rotamers <- length(names)
  
  #binary matrix indicating if each amino acid is present or not at position i in the sequence j
  AA<-lapply(1:nb_angles, function(i){
    t(table(c(rotamers[,i], names), row.names=c(1:(nb_frames+nb_rotamers))))[1:nb_frames,]
  })
  
  COV2 <- matrix(0, ncol= nb_angles, nrow= nb_angles)
  
  #Setting columns and rows names before matrix reduction
  rownames(COV2) <- pos_names
  colnames(COV2) <- pos_names
  
  # Calculating MI score for each valid position
  for(i in 1:nb_angles){
    mat_i <- AA[[i]] #matrix nb_frames*nb_rotamers
    
    cat(paste("pos_i : ", i, "\n"))
    
    for(j in i:nb_angles){
      mat_j <- AA[[j]] #matrix nb_frames*nb_rotamers
      
      #sequences in the alignment without gapped residues at positions i and k (no "-" in the positions i and k)
      Valid <- which(AA[[i]]%*%matrix(1, nrow= nb_rotamers)*AA[[j]]%*%matrix(1, nrow= nb_rotamers)!=0)
      
      #the number of sequences in the alignment without gapped residues at positions i and j (no "-" in the positions i and j)
      n <- length(Valid)
      
      Ex <- matrix(as.vector(t(AA[[i]])%*%matrix(1,nrow= nb_frames)), ncol= nb_rotamers, nrow= nb_rotamers)*t(matrix(as.vector(t(AA[[j]])%*%matrix(1,nrow= nb_frames)), ncol= nb_rotamers, nrow= nb_rotamers))/n
      Obs <- t(AA[[i]])%*%AA[[j]]
      COV2[i, j] <- sum((Obs-Ex)*(Obs-Ex)/n)
    }
  }
  
  COV2<-COV2+t(COV2) #Complete the second triangular part of the matrix
  diag(COV2) <- 0

  #Special square calculation for rotamers
  COV2 <- COV2^2
  
  #Removing "Inf" and "NA" values
  COV2[is.infinite(COV2)] <- 0
  COV2[is.na(COV2)] <- 0
  
  selected_pos_names <- colnames(COV2)
  if(!is.null(res_selection)){
    #Position selection
    tor.seq <- dynamic_struct$tor.seq
    
    residue_selection.inds <- which(tor.seq %in% res_selection)
    COV2 <- COV2[, residue_selection.inds]
    COV2 <- COV2[residue_selection.inds, ]
    
    nb_angles <- length(COV2[1,])
    selected_pos_names <- colnames(COV2)
  }
  
  #Removing autocorrelations
  angles_num <- c()
  if(auto_pairing == FALSE){
    ballesteros <- grepl("chi", selected_pos_names, fixed= TRUE)
    ballesteros <- unlist(lapply(1:length(selected_pos_names), function(i){if(ballesteros[i] == TRUE) 1 else 0}))
    nb_ballesteros <- sum(ballesteros)
    if(nb_ballesteros != 0){
      angles_num <- unlist(lapply(selected_pos_names, function(x){strsplit(x, "[.]")[[1]][1]}))
    } else {
      angles_num <- selected_pos_names
    }
    
    for(i in 1:nb_angles){
      for(j in i:nb_angles){
	ang_i <- angles_num[i]
	ang_j <- angles_num[j]
	
	if(ang_i == ang_j){
	  COV2[i,j] <- 0
	  COV2[j,i] <- 0
	}
      }
    }
  }
  
  res <- list()
  res$gross <- COV2
  
  if(z_score){
    COV2 <- (COV2-mean(COV2))/sd(COV2)
  }
  
  res$normalized <- COV2
	
  return(res)
}
