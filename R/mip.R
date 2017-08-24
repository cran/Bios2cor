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

mip <- function (align, fileHelix= NULL, diag= 0, fileCSV = NULL, gap_val= 0.8, z_score= TRUE) {

	if(!is.null(fileHelix)){
		#each line of the helix file looks like : "name1:x1-x2-x3"
		Helix<-read.table(fileHelix, sep= "-", stringsAsFactors= FALSE)
		
		#so we split each line like this : "name1:x1" "x2" "x3" (read.table can only have one separator)
		#then we split the "name1:x1" : "name1" "x1" for each line of our data frame Helix we keep "x1"
		lapply(1:length(Helix[,1]), function(x){ print(Helix[x,1]); Helix[x,1] <<- strsplit(Helix[x,1], split= ":")[[1]][2] })
	}
	
	msa<-align
	N <- matrix(as.vector(unlist(msa)), ncol = length(msa[[1]]),byrow = TRUE)
	nb_pos <- length(N[1,]) #number of positions in the alignment
	nb_seq <- length(N[,1]) #number of sequences in the alignment

	if(!is.null(fileHelix)){
		helix_inds <- unlist(lapply(1:7,function(i){Helix[i,1]:Helix[i,3]}))
		no_helix_inds <- setdiff(1:nb_pos, helix_inds)
		
		helix <- unlist(lapply(1:7,function(i){paste(rep(i,length(Helix[i,1]:Helix[i,3])),(Helix[i,1]:Helix[i,3]+50-Helix[i,2]),sep=".")}))
		
# 		colnames(N) <- rep(0.0, nb_pos)
		colnames(N) <- c(1:nb_pos)
		colnames(N)[helix_inds] <- helix
		
# 		N <- N[,unlist(lapply(1:7,function(i){Helix[i,1]:Helix[i,3]}))]
# 		nb_pos <- length(N[1,]) #number of positions in the alignment
# 		nb_seq <- length(N[,1]) #number of sequences in the alignment
# 	
# 		colnames(N)<-unlist(lapply(1:7,function(i){paste(rep(i,length(Helix[i,1]:Helix[i,3])),(Helix[i,1]:Helix[i,3]+50-Helix[i,2]),sep=".")}))
	}
	else {
		colnames(N)<-c(1:nb_pos)
	}
	names<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","-")
	AA<-lapply(1:length(N[1,]),function(i){t(table(c(N[,i],names),row.names=c(1:(length(N[,i])+21))))[1:length(N[,i]),-1]})

	names<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
	nb_aa <- length(names)
	
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
	
	MI <- matrix(0, ncol= nb_pos, nrow= nb_pos)
	
	freq_ij <- matrix(0, ncol= nb_aa, nrow= nb_aa, dimnames=list(names, names))
	
	# Calculating MI score for each valid position
	for(i in 1:nb_valid_pos){
	    pos_i <- Valid_pos[i] #current valid position
	    mat_i <- AA[[pos_i]] #matrix nb_seq*nb_aa
	    
	    cat(paste("pos_i : ", pos_i, "\n"))
	    
	    #amino acids frequency at position i in the alignment
	    freq_i <- colSums(mat_i)/nb_seq
	    
	    for(j in i:nb_valid_pos){
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
	
	#Correction P
	P <- matrix(0, ncol= nb_pos, nrow= nb_pos)
	
	mean_cov=sum(MI)/(length(which(MI!=0)))
	for(i in 1:(nb_valid_pos-1)){
		pos_i <- Valid_pos[i] #current valid position
		mean_i<-(sum(MI[pos_i,])+sum(MI[,pos_i]))/(nb_pos-1)

		for(k in (i+1):nb_valid_pos){
			pos_k <- Valid_pos[k] #current valid position
			mean_k<-(sum(MI[pos_k,])+sum(MI[,pos_k]))/(nb_pos-1)
			P[pos_i, pos_k]<-((mean_i*mean_k))/(mean_cov)
		}
	}

	diag(P)<-0

	COV2<-matrix(0, ncol= nb_pos, nrow= nb_pos)
	COV2 <- MI-P #MIP final value
	
	#Setting columns and rows names before matrix reduction
	rownames(COV2)<-paste(rep("pos", nb_pos), colnames(N),sep="")
	colnames(COV2)<-paste(rep("pos", nb_pos), colnames(N),sep="")
	
	COV2 <- COV2+t(COV2) #Complete the second triangular part of the matrix
	diag(COV2) <- diag
	
	#Reducting the final correlation matrix to the valid positions
	COV2 <- COV2[Valid_pos,]
	COV2 <- COV2[,Valid_pos]
	
	#Removing "Inf" and "NaN" values
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

