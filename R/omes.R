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
omes <- function(align, fileHelix = NULL , diag= 0, fileCSV = NULL, gap_val= 0.8, z_score= TRUE) {

    if(!is.null(fileHelix)){
	#each line of the helix file looks like : "name1:x1-x2-x3"
	Helix<-read.table(fileHelix, sep= "-", stringsAsFactors= FALSE)
	
	#so we split each line like this : "name1:x1" "x2" "x3" (read.table can only have one separator)
	#then we split the "name1:x1" : "name1" "x1" for each line of our data frame Helix we keep "x1"
	lapply(1:length(Helix[,1]), function(x){ print(Helix[x,1]); Helix[x,1] <<- strsplit(Helix[x,1], split= ":")[[1]][2] })
    }
    ams<-align
    N <- matrix(as.vector(unlist(ams)), ncol= length(ams[[1]]), byrow= TRUE)
    nb_pos <- length(N[1,]) #number of positions in the alignment
    nb_seq <- length(N[,1]) #number of sequences in the alignment
    
    if(!is.null(fileHelix)){
      helix_inds <- unlist(lapply(1:7,function(i){Helix[i,1]:Helix[i,3]}))
      no_helix_inds <- setdiff(1:nb_pos, helix_inds)
      
      helix <- unlist(lapply(1:7,function(i){paste(rep(i,length(Helix[i,1]:Helix[i,3])),(Helix[i,1]:Helix[i,3]+50-Helix[i,2]),sep=".")}))
      
      colnames(N) <- c(1:nb_pos)
      colnames(N)[helix_inds] <- helix
    }
    else {
      colnames(N)<-c(1:nb_pos)
    }
    names<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","-")

    #binary matrix indicating if each amino acid is present or not at position i in the sequence j
    AA<-lapply(1:nb_pos, function(i){
	    t(table(c(N[,i],names),row.names=c(1:(nb_seq+21))))[1:nb_seq, -1]
    })
    
    names<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
    nb_aa <- length(names)
    
    COV2<-matrix(0, ncol= nb_pos, nrow= nb_pos)
    
    #Setting columns and rows names before matrix reduction
    rownames(COV2)<-paste(rep("pos", nb_pos), colnames(N), sep="")
    colnames(COV2)<-paste(rep("pos", nb_pos), colnames(N), sep="")
    
    #Contains sequences that have unvalid gap proportion (> gap_val)
    Unvalid_pos <- c()
    
    #Looking for positions that have a correct gap proportion (according to "gap_val" argument)
    for(current_pos in 1:nb_pos){
      nb_gap_i <- nb_seq-sum(AA[[current_pos]]) #number of gaps in position "current_pos"
      ratio_i <- nb_gap_i/nb_seq #proportion of gaps in position "current_pos"
      if(ratio_i > gap_val) Unvalid_pos <- c(Unvalid_pos, current_pos)
    }
    nb_unvalid_pos <- length(Unvalid_pos)
    
    Valid_pos <- setdiff(1:nb_pos, Unvalid_pos)
    nb_valid_pos <- length(Valid_pos)
    
    # Calculating MI score for each valid position
    for(i in 1:nb_valid_pos){
      pos_i <- Valid_pos[i] #current valid position
      mat_i <- AA[[pos_i]] #matrix nb_seq*nb_aa
      
      cat(paste("pos_i : ", pos_i, "\n"))
      
      for(j in i:nb_valid_pos){
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
    
    COV2<-COV2+t(COV2) #Complete the second triangular part of the matrix
    diag(COV2) <- diag
    
    #Reducting the final correlation matrix to the valid positions
    COV2 <- COV2[Valid_pos,]
    COV2 <- COV2[,Valid_pos]

    #Removing "Inf" values
    COV2[is.infinite(COV2)] <- 0
    COV2[is.na(COV2)] <- 0
    
    res <- list()
    res$gross <- COV2
    
    if(z_score){
      COV2 <- (COV2-mean(COV2))/sd(COV2)
    }
    
    res$normalized <- COV2
    
    
    return(res)
}


