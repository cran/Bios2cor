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

dynamic_omes <- function(dynamic_structure, rotamers, res_selection= c("C","I","L","M","V","R","H","K","D","E","N","Q","F","Y","W","T","S","P")){
  
  if(missing(dynamic_structure)){
    stop("An object of class 'structure' is required as first argument")
  }

  if (missing(rotamers)){
    stop("A matrix of type 'rotamers' is required as second argument")
  }

  nb_angles <- length(rotamers[1,])
  nb_frames <- length(rotamers[,1])
  
  pos_names <- colnames(rotamers)
  
  names <- unique(as.vector(rotamers))
  nb_rotamers <- length(names)
  
 # Binary matrix indicating which rotamer is present or not at dihedral position i for each trajectory frame 
  ROT<-lapply(1:nb_angles, function(i){
    t(table(c(rotamers[,i], names), row.names=c(1:(nb_frames+nb_rotamers))))[1:nb_frames,]
  })
  
  COV2 <- matrix(0, ncol= nb_angles, nrow= nb_angles)
  
  #Setting columns and rows names before matrix reduction
  rownames(COV2) <- pos_names
  colnames(COV2) <- pos_names
  
  # Calculating OMES score for each position
  for(i in 1:nb_angles){
    mat_i <- ROT[[i]] #matrix nb_frames*nb_rotamers
    
    cat(paste("pos_i : ", i, "\n"))
    
    for(j in i:nb_angles){
      mat_j <- ROT[[j]] #matrix nb_frames*nb_rotamers   
      Valid <- which(ROT[[i]]%*%matrix(1, nrow= nb_rotamers)*ROT[[j]]%*%matrix(1, nrow= nb_rotamers)!=0)    #possible pairs
      n <- length(Valid)
      Ex <- matrix(as.vector(t(ROT[[i]])%*%matrix(1,nrow= nb_frames)), ncol= nb_rotamers, nrow= nb_rotamers)*t(matrix(as.vector(t(ROT[[j]])%*%matrix(1,nrow= nb_frames)), ncol= nb_rotamers, nrow= nb_rotamers))/n
      Obs <- t(ROT[[i]])%*%ROT[[j]]
      COV2[i, j] <- sum((Obs-Ex)*(Obs-Ex)/n)
    }
  }
  
  COV2<-COV2+t(COV2) #Complete the second triangular part of the matrix
  diag(COV2) <- 0

  #Removing "Inf" and "NA" values
  COV2[is.infinite(COV2)] <- 0
  COV2[is.na(COV2)] <- 0
  
  selected_pos_names <- colnames(COV2)

  if(!is.null(res_selection)){
    #Position selection
    tor.seq <- dynamic_structure$tor.seq
    
    residue_selection.inds <- which(tor.seq %in% res_selection)
    COV2 <- COV2[, residue_selection.inds]
    COV2 <- COV2[residue_selection.inds, ]
    
    nb_angles <- length(COV2[1,])
    selected_pos_names <- colnames(COV2)
  }
  

 
res <- list()
    # save matrix of scores 
    res$score <- COV2
    
    # Compute and save matrix of Z_scores
    # Mean and stdev must be calculated on off diagonal elements  
    mean_up <- mean(COV2[upper.tri(COV2)])
    stdev <- sd(COV2[upper.tri(COV2)])
    COV3 <- (COV2-mean_up)/stdev
    diag(COV3) <- 0
    res$Zscore <- COV3

    # Save matrices of scores and Zscores without auto correlation
    # Create and save matrix with 0 values for autocorrelation (dihedral angles within the same residue)
    COV5 <- COV2
    res_num <- c()
    res_num <- unlist(lapply(selected_pos_names, function(x){strsplit(x, "[.]")[[1]][1]}))

    # Create matrix with NA that will be used for Zscores 
    for(i in 1:nb_angles){
      for(j in i:nb_angles){
	res_i <- res_num[i]
	res_j <- res_num[j]
	
	if(res_i == res_j){
	  COV5[i,j] <- NA
	  COV5[j,i] <- NA
	}
      }
    }
 


    # Save matrix of score with 0 for autocorrelation  
    COV4 <-COV5
    COV4[is.na(COV4)] <- 0
    res$score_noauto <- COV4

    # Compute and save matrix of Z_scores without auto correlation  
    # Mean and stdev are calculated on non NA elements of upper triangle 
    mean_noauto <- mean(COV5[upper.tri(COV5)],na.rm=TRUE)
    stdev_noauto <- sd(COV5[upper.tri(COV5)],na.rm=TRUE)
    
    COV5[is.na(COV5)] <- 0
    COV5 <- (COV5-mean_noauto)/stdev_noauto

    # Correct diag and autocorrelated elements	
    diag(COV5) <- 0
    
    for(i in 1:nb_angles){
      for(j in i:nb_angles){
	res_i <- res_num[i]
	res_j <- res_num[j]
	
	if(res_i == res_j){
	  COV5[i,j] <- 0
	  COV5[j,i] <- 0
	}
      }
    }
  
    res$Zscore_noauto <- COV5

  return(res)
}
