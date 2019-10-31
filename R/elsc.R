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
elsc <- function (align, gap_ratio = 0.2){

  if ((gap_ratio < 0) | (gap_ratio > 1)) {
    stop("Error in elsc argument: gap_ratio must be in the [0,1] range")
  }

    diag <- 0

    msa<-align
    MSA <- matrix(as.vector(unlist(msa)), ncol = length(msa[[1]]),byrow = TRUE)
    nb_pos <- length(MSA[1,]) #number of positions in the alignment
    nb_seq <- length(MSA[,1]) #number of sequences in the alignment
    colnames(MSA)<-c(1:nb_pos)
    pos_names <- colnames(MSA)  

    gap <- 1-gap_ratio      #gap value indicates the minimal ratio of aa to nb_seq in the MSA 
    if (gap < 1/nb_seq) {
      gap <- 1/nb_seq     # positions must have at leat ONE aa to be taken into account (removes gap column)
    }

    names<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","-")

    AA<-lapply(1:nb_pos,function(i){
	t(table(c(MSA[,i],names),row.names=c(1:(nb_seq+21))))[1:nb_seq, -1]
    })

    names<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
    nb_aa <- length(names)
	
    COV2<-matrix(0, ncol= nb_pos, nrow= nb_pos)
	
    # Setting columns and rows names before matrix reduction
    rownames(COV2)<-pos_names
    colnames(COV2)<-pos_names
	
    # Determining valid positions with gap ratio under the (1 - gap) limit 
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
	
    # Calculating ELSC score for each valid position
    for(i in 1:nb_Valid_pos){
        pos_i <- Valid_pos[i] # current valid position
        cat(paste("pos_i : ", pos_i, "\n"))
	for(j in i:nb_Valid_pos){
	    pos_j <- Valid_pos[j]
		
	    aln_tot_i<-AA[[pos_i]]
	    aln_tot_i[is.na(aln_tot_i)] <- 0
		
	    aln_tot_j<-AA[[pos_j]]
	    aln_tot_j[is.na(aln_tot_j)] <- 0
		
	    Nyj<-colSums(aln_tot_j) # number of amino acids y in position j in the full alignment
	    v<-sort(colSums(aln_tot_i),decreasing = TRUE)
	    v<-unique(c(which(aln_tot_i[,names(v[1])] == 1))) # sequences where the most occuring amino acid in position i appears 

	    ss_aln_i<-aln_tot_i[v,]
	    nb_seq_ss_aln_i <- length(ss_aln_i[,1])
		
	    ss_aln_j<-aln_tot_j[v,]
	    nb_seq_ss_aln_j <- length(ss_aln_j[,1])
		
	    n<-colSums(ss_aln_j) # number of amino acids y in position j in the sub-alignment
	    mf<-(Nyj/nb_seq)*nb_seq_ss_aln_j
		
	    calc<-matrix(0, ncol=2, nrow= nb_aa, dimnames = list(names,c("r", "m")))
	    calc[,2]<-trunc(mf) # inferior mf round
	    calc[,1]<-(mf-calc[,2]) # remainder

	    calc<-calc[order(calc[,1],rownames(calc), decreasing =T),] # order m by decreasing order of r then by alpha inverse order of aa
		
	    # The sum of the "m" column must be equal to the sum of "n"
	    if(sum(calc[,2]) > sum(n)){
		z<-nb_aa
		while((sum(calc[,2]))!=(sum(n))){ 
		    calc[z,2]<-(calc[z,2]-1)
		    if(z == 1) z <- nb_aa else z<-z-1
		}
	    } else {
		z<-1
		while((sum(calc[,2]))!=(sum(n))){
		    calc[z,2]<-(calc[z,2]+1)
		    if(z == nb_aa) z <- 1 else z<-z+1
		}
	    }
		
	    m<-calc[order(rownames(calc)),2]
	    COV2[pos_i, pos_j]<-(-log(prod(choose(Nyj, n)/choose(Nyj, m)))) #matrice asymetrique - vrais scores obtenus 
	    }
	}



	#Complete matrix second triangular part
	COV2 <- COV2 + t(COV2)
	diag(COV2) <- diag

	#Reducting the final correlation matrix to the valid positions
	COV2 <- COV2[Valid_pos,]
	COV2 <- COV2[,Valid_pos]

	#Removing "Inf" and "NaN" values
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
