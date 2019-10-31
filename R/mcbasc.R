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
mcbasc <- function (align, gap_ratio = 0.2) {

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




    matrice_mclachlan <- matrix(c(
			8,1,3,4,1,3,3,2,3,2,3,3,4,3,2,4,3,3,1,1,		
			1,9,1,0,0,1,3,1,0,0,3,1,0,0,1,2,2,1,2,1,		
			3,1,8,5,1,3,4,0,3,1,2,5,3,4,1,3,3,1,0,1,		
			4,0,5,8,0,3,2,1,4,1,1,4,4,5,3,4,4,2,1,2,		
			1,0,1,0,9,0,4,3,0,5,5,0,1,0,1,2,1,3,6,6,		
			3,1,3,3,0,8,2,1,3,1,1,3,3,2,3,3,2,2,1,0,		
			3,3,4,2,4,2,8,2,4,2,3,4,3,4,5,3,4,2,3,4,		
			2,1,0,1,3,1,2,8,1,5,5,1,1,0,1,2,3,5,3,3,		
			3,0,3,4,0,3,4,1,8,2,1,4,3,4,5,3,3,2,1,1,		
			2,0,1,1,5,1,2,5,2,8,6,1,1,3,2,2,3,5,3,3,		
			3,3,2,1,5,1,3,5,1,6,8,2,1,3,1,2,3,4,1,2,		
			3,1,5,4,0,3,4,1,4,1,2,8,1,4,3,5,3,1,0,2,		
			4,0,3,4,1,3,3,1,3,1,1,1,8,3,3,3,3,2,0,0,		 
			3,0,4,5,0,2,4,0,4,3,3,4,3,8,5,4,3,2,2,1,		
			2,1,1,3,1,3,5,1,5,2,1,3,3,5,8,4,3,2,3,2,		
			4,2,3,4,2,3,3,2,3,2,2,5,3,4,4,8,5,2,3,3,		
			3,2,3,4,1,2,4,3,3,3,3,3,3,3,3,5,8,3,2,1,		
			3,1,1,2,3,2,2,5,2,5,4,1,2,2,2,2,3,8,2,3,		
			1,2,0,1,6,1,3,3,1,3,1,0,0,2,3,3,2,2,9,6,		
			1,1,1,2,6,0,4,3,1,3,2,2,0,1,2,3,1,3,6,9),nrow = 20)
	
    names<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
    rownames(matrice_mclachlan) <- names
    colnames(matrice_mclachlan) <- names

    A<-lapply(1:nb_pos, function(i){
	  matrix(0, ncol=nb_seq, nrow=nb_seq, dimnames = list(c(MSA[,i]), c(MSA[,i])))
    })

    for (k in 1:nb_pos) {
 	mat<-A[[k]]
	for (i in names) {
		for (j in names) {
			mat[rownames(mat)==i, colnames(mat)==j] <- matrice_mclachlan[rownames(matrice_mclachlan)==i, colnames(matrice_mclachlan)==j]
		}
	}
	A[[k]] <- mat
    }

    COV2<-matrix(0, ncol= nb_pos, nrow= nb_pos)
	
    # Setting columns and rows names before matrix reduction
    rownames(COV2)<-pos_names
    colnames(COV2)<-pos_names
	
    # Calculating McBASC score for each valid position
    for(i in 1:nb_Valid_pos){
	pos_i <- Valid_pos[i] #current valid position
	mat_i <- A[[pos_i]] #matrix nb_seq*nb_aa
	    
	cat(paste("pos_i : ", pos_i, "\n"))
	    
 	for(j in i:nb_Valid_pos){
	    pos_j <- Valid_pos[j]
	    mat_j <- A[[pos_j]] #matrix nb_seq*nb_aa
		
	    SUM<-sum((mat_i-mean(mat_i))*(mat_j-mean(mat_j)))
	    COV2[pos_i, pos_j] <- (1/((nb_seq)^2))*(SUM/(sd(mat_i)*sd(mat_j)))
	}
    }


    # Complete the second triangular part of the matrix	
    COV2 <- COV2+t(COV2)
    diag(COV2) <- diag
	
    # Reducting the final correlation matrix to the valid positions
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

