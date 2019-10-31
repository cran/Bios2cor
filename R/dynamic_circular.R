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

dynamic_circular <- function(dynamic_structure, res_selection= c("C","I","L","M","V","R","H","K","D","E","N","Q","F","Y","W","T","S","P")){

  if(missing(dynamic_structure)){
    stop("A 'dynamic_structure' object is required")
  }

  #  Importing torsional information
  tor <- dynamic_structure$tor
  nb_torsions <- length(tor[1,])
  
  # Torsional angles
  sidechain.tor <- tor[, 1:nb_torsions]
  
  dihed_names <- colnames(sidechain.tor)
  
  # Transforms the torsion matrix in circular object (circular package)
  dihed <- circular(sidechain.tor, units="degrees", type="angles", modulo="2pi")

  # Computes circular correlation
  dihed_corr <- cor.circular(dihed)

  colnames(dihed_corr) <- dihed_names
  rownames(dihed_corr) <- dihed_names
  
  diag(dihed_corr) <- 0
  
  # Squares because only absolute value of correlation matters
  dihed_corr <- dihed_corr^2    
  
  dihed_corr_names <- colnames(dihed_corr)
  
  if(!is.null(res_selection)){
    # Position selection
    tor.seq <- dynamic_structure$tor.seq
    
    residue_selection.inds <- which(tor.seq %in% res_selection)
    dihed_corr <- dihed_corr[, residue_selection.inds]
    dihed_corr <- dihed_corr[residue_selection.inds, ]
    
    nb_angles  <- length(dihed_corr[1,])
    selected_dihed_names <- colnames(dihed_corr)
  }
  

  COV2 <- dihed_corr

  # Removing "Inf" and "NA" values
  COV2[is.infinite(COV2)] <- 0
  COV2[is.na(COV2)] <- 0


res <- list()
    # Save matrix of scores 
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
    res_num <- unlist(lapply(selected_dihed_names, function(x){strsplit(x, "[.]")[[1]][1]}))

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

    ### Correct diag and autocorrelated elements	
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
