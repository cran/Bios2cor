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
omes <- function(align, gap_ratio = 0.2) {

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

    # Binary matrix indicating which amino acid is present or not at position i in the sequence j
    AA<-lapply(1:nb_pos, function(i){
	    t(table(c(MSA[,i],names),row.names=c(1:(nb_seq+21))))[1:nb_seq, -1]
    })
    
    names<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
    nb_aa <- length(names)
    
    COV2<-matrix(0, ncol= nb_pos, nrow= nb_pos)
    
    # Setting columns and rows names before matrix reduction
    rownames(COV2)<-pos_names
    colnames(COV2)<-pos_names 
    
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
    
    # Calculating score for each valid position
    for(i in 1:nb_Valid_pos){
      pos_i <- Valid_pos[i] #current valid position
      mat_i <- AA[[pos_i]] #matrix nb_seq*nb_aa
      
      cat(paste("pos_i : ", pos_i, "\n"))
      
      for(j in i:nb_Valid_pos){
	pos_j <- Valid_pos[j]
	mat_j <- AA[[pos_j]] #matrix nb_seq*nb_aa
	
	#sequences in the alignment without gapped residues at positions i and k (no "-" in the positions i and k)
	Valid<-which(AA[[pos_i]]%*%matrix(1, nrow= nb_aa)*AA[[pos_j]]%*%matrix(1, nrow= nb_aa)!=0)
	
	#the number of sequences in the alignment without gapped residues at positions i and j (no "-" in the positions i and j)
	n<-length(Valid)
	
	Ex<-matrix(as.vector(t(AA[[pos_i]])%*%matrix(1,nrow= nb_seq)), ncol= nb_aa, nrow= nb_aa)*t(matrix(as.vector(t(AA[[pos_j]])%*%matrix(1,nrow= nb_seq)), ncol= nb_aa, nrow= nb_aa))/n
	Obs<-t(AA[[pos_i]])%*%AA[[pos_j]]
	COV2[pos_i, pos_j]<-sum((Obs-Ex)*(Obs-Ex)/n)
      }
    }
    
    #Complete the second triangular part of the matrix
    COV2<-COV2+t(COV2) 
    diag(COV2) <- diag
    
    #Reducting the final correlation matrix to the valid positions
    COV2 <- COV2[Valid_pos,]
    COV2 <- COV2[,Valid_pos]

    #Removing "Inf" values
    COV2[is.infinite(COV2)] <- 0
    COV2[is.na(COV2)] <- 0
    
    res <- list()
  
    # save matrix of score 
    res$score <- COV2
    
    # compute and save matrix of Z_scores
    # mean and stdev must be calculated on off diagonal elements  
    mean_up <- mean(COV2[upper.tri(COV2)])
    stdev <- sd(COV2[upper.tri(COV2)])
    
    COV2 <- (COV2-mean_up)/stdev
    diag(COV2) <- diag
       
    res$Zscore <- COV2
    
    return(res)
}


