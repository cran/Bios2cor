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
rotamer_mip <- function(dynamic_struct, rotamers, res_selection= c("C","I","L","M","V","R","H","K","D","E","N","Q","F","Y","W","T","S","P"), z_score= TRUE, auto_pairing= FALSE){
  nb_angles <- length(rotamers[1,])
  nb_frames <- length(rotamers[,1])
  
  pos_names <- colnames(rotamers)
  
  names <- unique(as.vector(rotamers))
  nb_rotamers <- length(names)
  
  #binary matrix indicating if each amino acid is present or not at position i in the sequence j
  AA<-lapply(1:nb_angles, function(i){
      t(table(c(rotamers[,i], names), row.names=c(1:(nb_frames+nb_rotamers))))[1:nb_frames,]
  })
  
  MI <- matrix(0, ncol= nb_angles, nrow= nb_angles)
	
  freq_ij <- matrix(0, ncol= nb_rotamers, nrow= nb_rotamers, dimnames=list(names, names))
  
  # Calculating MI score for each valid position
  for(i in 1:nb_angles){
      mat_i <- AA[[i]] #matrix nb_frames*nb_rotamers
      
      cat(paste("i : ", i, "\n"))
      
      #amino acids frequency at position i in the alignment
      freq_i <- colSums(mat_i)/nb_frames
      
      for(j in i:nb_angles){
	  mat_j <- AA[[j]] #matrix nb_seq*nb_rotamers
	  
	  #amino acids frequency at position j in the alignment
	  freq_j <- colSums(mat_j)/nb_frames #matrix 1*nb_rotamers

	  #amino acids frequency at positions i AND j in the alignment
	  freq_p <- ((t(mat_i))%*%mat_j)/nb_frames #matrix nb_rotamers*nb_rotamers
	  
	  for (k in names) {
	      for (l in names) {
		  freq_ij[k,l]<-freq_i[k]*freq_j[l]
	      }
	  }
	  
	  LOG<-(log(freq_p, base=400)-log(freq_ij, base=400))
	  LOG<-replace(LOG, which(LOG=="-Inf"),0)
	  LOG<-replace(LOG, which(LOG=="Inf"),0)
	  LOG<-replace(LOG, which(LOG=="NaN"),0)
	  
	  MI[i, j]<-sum((freq_p)*(LOG))
      }
  }
  
  #Correction P
  P <- matrix(0, ncol= nb_angles, nrow= nb_angles)
  
  mean_cov=sum(MI)/(length(which(MI!=0)))
  for(i in 1:(nb_angles-1)){
	  mean_i<-(sum(MI[i,])+sum(MI[,i]))/(nb_angles-1)

	  for(k in (i+1):nb_angles){
		  mean_k<-(sum(MI[k,])+sum(MI[,k]))/(nb_angles-1)
		  P[i, k]<-((mean_i*mean_k))/(mean_cov)
	  }
  }

  diag(P)<-0
  
  COV2<-matrix(0, ncol= nb_angles, nrow= nb_angles)
  COV2 <- MI-P #MIP final value
  
  #Setting columns and rows names before matrix reduction
  rownames(COV2)<-paste(pos_names, sep="")
  colnames(COV2)<-paste(pos_names, sep="")
  
  COV2<-COV2+t(COV2) #Complete the second triangular part of the matrix
  diag(COV2) <- 0
  
  #Special square calculation for rotamers
  COV2 <- COV2^2
    
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
 
