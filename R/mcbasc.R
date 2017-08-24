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
mcbasc <- function (align, fileHelix= NULL, diag= 0, fileCSV= NULL, gap_val= 0.8, z_score= TRUE) {

	if(!is.null(fileHelix)){
		#each line of the helix file looks like : "name1:x1-x2-x3"
		Helix<-read.table(fileHelix, sep= "-", stringsAsFactors= FALSE)
		
		#so we split each line like this : "name1:x1" "x2" "x3" (read.table can only have one separator)
		#then we split the "name1:x1" : "name1" "x1" for each line of our data frame Helix we keep "x1"
		lapply(1:length(Helix[,1]), function(x){ print(Helix[x,1]); Helix[x,1] <<- strsplit(Helix[x,1], split= ":")[[1]][2] })
	}
	
	ams<-align
	N <- matrix(as.vector(unlist(ams)), ncol = length(ams[[1]]),byrow = TRUE)
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
	  matrix(0, ncol=nb_seq, nrow=nb_seq, dimnames = list(c(N[,i]), c(N[,i])))
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
	
	#Setting columns and rows names before matrix reduction
	rownames(COV2)<-paste(rep("pos", nb_pos), colnames(N),sep="")
	colnames(COV2)<-paste(rep("pos", nb_pos), colnames(N),sep="")
	
	#Contains sequences that have unvalid gap proportion (> gap_val)
	Unvalid_pos <- c()
	
	#Looking for positions that have a correct gap proportion (according to "gap_val" argument)
	for(current_pos in 1:nb_pos){
	    nb_gap_i <- sum(rownames(A[[current_pos]]) == "-") #number of gaps in position "current_pos"
	    ratio_i <- nb_gap_i/nb_seq #proportion of gaps in position "current_pos"
	    if(ratio_i > gap_val) Unvalid_pos <- c(Unvalid_pos, current_pos)
	}
	nb_unvalid_pos <- length(Unvalid_pos)
	
	Valid_pos <- setdiff(1:nb_pos, Unvalid_pos)
	nb_valid_pos <- length(Valid_pos)
	
	# Calculating McBASC score for each valid position
	for(i in 1:nb_valid_pos){
	    pos_i <- Valid_pos[i] #current valid position
	    mat_i <- A[[pos_i]] #matrix nb_seq*nb_aa
	    
	    cat(paste("pos_i : ", pos_i, "\n"))
	    
	    for(j in i:nb_valid_pos){
		pos_j <- Valid_pos[j]
		mat_j <- A[[pos_j]] #matrix nb_seq*nb_aa
		
		somme<-sum((mat_i-mean(mat_i))*(mat_j-mean(mat_j)))
		COV2[pos_i, pos_j] <- (1/((nb_seq)^2))*(somme/(sd(mat_i)*sd(mat_j)))
	    }
	}
	
	COV2 <- COV2+t(COV2)
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

