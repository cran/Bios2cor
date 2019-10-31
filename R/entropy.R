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

entropy <- function (align, gap_ratio = 0.2) {

    if ((gap_ratio < 0) | (gap_ratio > 1)) {
        stop("Error in entropy argument: gap_ratio must be in the [0,1] range.")
    }

    MSA <- matrix(as.vector(unlist(align)), ncol= length(align[[1]]), byrow= TRUE)
    nb_pos <- length(MSA[1,]) #number of positions in the alignment
    nb_seq <- length(MSA[,1]) #number of sequences in the alignment
    colnames(MSA)<-c(1:nb_pos)

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

    entropy <- matrix(0, ncol= nb_pos, nrow= 1)

     # Calculating validity and entropy for each  position
    entropy <- unlist(lapply(1:nb_pos, function(i){
       mat_i <- AA[[i]] #matrix nb_seq*nb_aa
        S_i <- colSums(AA[[i]])
        Tot_i <- sum(S_i)
        if (Tot_i/nb_seq >= gap) {
      		percent_i <- S_i/Tot_i
      		entropy_i <- 0
      		for(j in 1:nb_aa){
         		if (S_i[[j]] == 0) {
				entropy_i <- entropy_i
	 		} else {
				entropy_i <- entropy_i - percent_i[[j]]*log(percent_i[[j]], base = 20) 
         		}
      		}	
      		entropy[i] <- entropy_i
    	} else {
		entropy[i] <- ''	
    	}
   })) 	
          
	names(entropy) <- c(1:nb_pos)
	return(entropy)
 

}
