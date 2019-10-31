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
dynamic_mi <- function(dynamic_structure, rotamers, res_selection= c("C","I","L","M","V","R","H","K","D","E","N","Q","F","Y","W","T","S","P")){

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
  
  MI <- matrix(0, ncol= nb_angles, nrow= nb_angles)
	
  freq_ij <- matrix(0, ncol= nb_rotamers, nrow= nb_rotamers, dimnames=list(names, names))
  
  # Calculating MI score for each position
  for(i in 1:nb_angles){
      mat_i <- ROT[[i]] # matrix nb_frames*nb_rotamers
      
      cat(paste("i : ", i, "\n"))
      
      # Rotamer frequency at dihedral position i in the trajectory
      freq_i <- colSums(mat_i)/nb_frames
      
      for(j in i:nb_angles){
	  mat_j <- ROT[[j]] #matrix nb_seq*nb_rotamers
	  
	  # Rotamer frequency at dihedral position j in the trajectory
	  freq_j <- colSums(mat_j)/nb_frames #matrix 1*nb_rotamers

	  # Rotamer frequency at dihedral positions i AND j in the trajectory
	  freq_p <- ((t(mat_i))%*%mat_j)/nb_frames #matrix nb_rotamers*nb_rotamers
	  
	  for (k in names) {
	      for (l in names) {
		  freq_ij[k,l]<-freq_i[k]*freq_j[l]
	      }
	  }
	  
	  LOG<-(log(freq_p, base=9)-log(freq_ij, base=9))
	  LOG<-replace(LOG, which(LOG=="-Inf"),0)
	  LOG<-replace(LOG, which(LOG=="Inf"),0)
	  LOG<-replace(LOG, which(LOG=="NaN"),0)
	  
	  MI[i, j]<-sum((freq_p)*(LOG))
      }
  }
  

  # MI final value
  COV2<-matrix(0, ncol= nb_angles, nrow= nb_angles)
  COV2 <- MI 
  
  # Setting columns and rows names before matrix reduction
  rownames(COV2)<-paste(pos_names, sep="")
  colnames(COV2)<-paste(pos_names, sep="")
  
  COV2<-COV2+t(COV2) #Complete the second triangular part of the matrix
  diag(COV2) <- 0
  
    
  # Removing "Inf" and "NA" values
  COV2[is.infinite(COV2)] <- 0
  COV2[is.na(COV2)] <- 0
 
  selected_pos_names <- colnames(COV2)
  
  if(!is.null(res_selection)){
    
    # Position selection
    tor.seq <- dynamic_structure$tor.seq
    
    residue_selection.inds <- which(tor.seq %in% res_selection)
    COV2 <- COV2[, residue_selection.inds]
    COV2 <- COV2[residue_selection.inds, ]
    
    nb_angles <- length(COV2[1,])
    selected_pos_names <- colnames(COV2)
  }
  
 
res <- list()
    # Save matrix of score 
    res$score <- COV2
    
    # Compute and save matrix of Z_scores
    # Mean and stdev must be calculated on off diagonal elements  
    mean_up <- mean(COV2[upper.tri(COV2)])
    stdev <- sd(COV2[upper.tri(COV2)])
    COV3 <- (COV2-mean_up)/stdev
    diag(COV3) <- 0
    res$Zscore <- COV3

    # Save matrices of score and Zscores without auto correlation
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
 
