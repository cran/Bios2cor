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

rotamer_circular <- function(dynamic_struct, res_selection= c("C","I","L","M","V","R","H","K","D","E","N","Q","F","Y","W","T","S","P"), z_score= TRUE, auto_pairing= FALSE){
  #Importing rotamer informations
  tor <- dynamic_struct$tor
  nb_residus <- length(tor[1,])
  
  #torsion angles
  lsidech.tor <- tor[, 1:nb_residus]
  
  pos_names <- colnames(lsidech.tor)
  
  #transforms the torsion matrix in circular object (circular package)
  dihed <- circular(lsidech.tor, units="degrees", type="angles", modulo="2pi")

  #computes circular correlation
  dihed_corr <- cor.circular(dihed)

  colnames(dihed_corr) <- pos_names
  rownames(dihed_corr) <- pos_names
  
  diag(dihed_corr) <- 0
  dihed_corr <- dihed_corr^2
  
  dihed_corr_names <- colnames(dihed_corr)
  
  if(!is.null(res_selection)){
    #Position selection
    tor.seq <- dynamic_struct$tor.seq
    
    residue_selection.inds <- which(tor.seq %in% res_selection)
    dihed_corr <- dihed_corr[, residue_selection.inds]
    dihed_corr <- dihed_corr[residue_selection.inds, ]
    
    nb_residus <- length(dihed_corr[1,])
    dihed_corr_names <- colnames(dihed_corr)
  }
  
  #Removing autocorrelations
  angles_num <- c()
  if(auto_pairing == FALSE){
    ballesteros <- grepl("chi", dihed_corr_names, fixed= TRUE)
    ballesteros <- unlist(lapply(1:length(dihed_corr_names), function(i){if(ballesteros[i] == TRUE) 1 else 0}))
    nb_ballesteros <- sum(ballesteros)
    if(nb_ballesteros != 0){
      angles_num <- unlist(lapply(dihed_corr_names, function(x){strsplit(x, "[.]")[[1]][1]}))
    } else {
      angles_num <- dihed_corr_names
    }
    
    for(i in 1:nb_residus){
      for(j in i:nb_residus){
	ang_i <- angles_num[i]
	ang_j <- angles_num[j]
	
	if(ang_i == ang_j){
	  dihed_corr[i,j] <- 0
	  dihed_corr[j,i] <- 0
	}
      }
    }
  }
  
  res <- list()
  res$gross <- dihed_corr

  if(z_score){
    dihed_corr <- (dihed_corr-mean(dihed_corr))/sd(dihed_corr)
  }
  
  res$normalized <- dihed_corr
  
  return(res)
}
