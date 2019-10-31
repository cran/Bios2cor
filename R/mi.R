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

mi <- function (align, gap_ratio = 0.2) {

    if ((gap_ratio < 0) | (gap_ratio > 1)) {
      stop("gap_ratio must be in the [0,1] range.")
    }

    diag <- 0	

    msa<-align
    MSA <- matrix(as.vector(unlist(msa)), ncol= length(msa[[1]]), byrow= TRUE)
    nb_pos <- length(MSA[1,]) #number of positions in the alignment
    nb_seq <- length(MSA[,1]) #number of sequences in the alignment
    colnames(MSA)<-c(1:nb_pos)
    pos_names <- colnames(MSA) 
	
    gap <- 1-gap_ratio      #gap value indicates the minimal ratio of aa to nb_seq in the MSA 
    if (gap < 1/nb_seq) {
      gap <- 1/nb_seq     # positions must have at leat ONE aa to be taken into account (removes gap column)
    }

    names<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","-")
    AA<-lapply(1:length(MSA[1,]),function(i){t(table(c(MSA[,i],names),row.names=c(1:(length(MSA[,i])+21))))[1:length(MSA[,i]),-1]})

    names<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
    nb_aa <- length(names)
	
    # Determining valid positions with correct gap ratio (equal to 1 - gap argument)
    Valid_pos <- c()
    for(i in 1:nb_pos){
       mat_i <- AA[[i]] #matrix nb_seq*nb_aa
        S_i <- colSums(AA[[i]])
        Tot_i <- sum(S_i)
        if (Tot_i/nb_seq >= gap) {
           Valid_pos <- c(Valid_pos, i)
        }
    }
    nb_Valid_pos <- length(Valid_pos)
	
    MI <- matrix(0, ncol= nb_pos, nrow= nb_pos)
    freq_ij <- matrix(0, ncol= nb_aa, nrow= nb_aa, dimnames=list(names, names))
	
    # Calculating MI score for each valid position
    for(i in 1:nb_Valid_pos){
	pos_i <- Valid_pos[i] #current valid position
	mat_i <- AA[[pos_i]] #matrix nb_seq*nb_aa
	    
	cat(paste("pos_i : ", pos_i, "\n"))
	    
	#amino acids frequency at position i in the alignment
	freq_i <- colSums(mat_i)/nb_seq
	    
	for(j in i:nb_Valid_pos){
	    pos_j <- Valid_pos[j]
	    mat_j <- AA[[pos_j]] #matrix nb_seq*nb_aa
		
	    #amino acids frequency at position j in the alignment
	    freq_j <- colSums(mat_j)/nb_seq #matrix 1*nb_aa

	    #amino acids frequency at positions i AND j in the alignment
	    freq_p <- ((t(mat_i))%*%mat_j)/nb_seq #matrix nb_aa*nb_aa
		
	    for (k in names) {
	        for (l in names) {
	    	    freq_ij[k,l]<-freq_i[k]*freq_j[l]
		}
	    }
		
	LOG<-(log(freq_p, base=400)-log(freq_ij, base=400))
	LOG<-replace(LOG, which(LOG=="-Inf"),0)
	LOG<-replace(LOG, which(LOG=="NaN"),0)
		
	MI[pos_i, pos_j]<-sum((freq_p)*(LOG))
	}
    }
	

    COV2<-matrix(0, ncol= nb_pos, nrow= nb_pos)
	
    # MI final value
    COV2 <- MI 
	
    # Setting columns and rows names before matrix reduction
    rownames(COV2)<-pos_names
    colnames(COV2)<-pos_names
	
    COV2 <- COV2+t(COV2) #Complete the second triangular part of the matrix
    diag(COV2) <- diag
	
    # Reduction of the final correlation matrix to the valid positions
    COV2 <- COV2[Valid_pos,]
    COV2 <- COV2[,Valid_pos]
	
    # Removing "Inf" and "NaN" values
    COV2[is.infinite(COV2)] <- 0
    COV2[is.na(COV2)] <- 0
	
    res <- list()
    res$score <- COV2
	
    # Compute and save matrix of Z_scores
    # Mean and stdev must be calculated on off diagonal elements  
    mean_up <- mean(COV2[upper.tri(COV2)])
    stdev <- sd(COV2[upper.tri(COV2)])
    
    COV2 <- (COV2-mean_up)/stdev
    diag(COV2) <- diag
       
    res$Zscore <- COV2
	
    return(res)

}

