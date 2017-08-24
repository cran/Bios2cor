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
elsc <- function (align, fileHelix= NULL, diag= 0, fileCSV= NULL, gap_val= 0.8, double_passing= FALSE, z_score= TRUE) {

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

	AA<-lapply(1:nb_pos,function(i){
		t(table(c(N[,i],names),row.names=c(1:(nb_seq+21))))[1:nb_seq, -1]
	})

	names<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
	nb_aa <- length(names)
	
	COV2<-matrix(0, ncol= nb_pos, nrow= nb_pos)
	
	#Setting columns and rows names before matrix reduction
	rownames(COV2)<-paste(rep("pos", nb_pos), colnames(N),sep="")
	colnames(COV2)<-paste(rep("pos", nb_pos), colnames(N),sep="")
	
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
	
	# Calculating ELSC score for each valid position
	for(i in 1:nb_valid_pos){
	    pos_i <- Valid_pos[i] #current valid position
	    
	    cat(paste("pos_i : ", pos_i, "\n"))
	    
	    for(j in i:nb_valid_pos){
		pos_j <- Valid_pos[j]
		
		aln_tot_i<-AA[[pos_i]]
		aln_tot_i[is.na(aln_tot_i)] <- 0
		
		aln_tot_j<-AA[[pos_j]]
		aln_tot_j[is.na(aln_tot_j)] <- 0
		
		Nyj<-colSums(aln_tot_j) #number of amino acid y in position j in the full alignment
		v<-sort(colSums(aln_tot_i),decreasing = TRUE)
		v<-unique(c(which(aln_tot_i[,names(v[1])] == 1))) #sequences where appears the most occuring amino acid in position i

		ss_aln_i<-aln_tot_i[v,]
		nb_seq_ss_aln_i <- length(ss_aln_i[,1])
		
		ss_aln_j<-aln_tot_j[v,]
		nb_seq_ss_aln_j <- length(ss_aln_j[,1])
		
		n<-colSums(ss_aln_j) #number of amino acid y in position j in the sub-alignment
		mf<-(Nyj/nb_seq)*nb_seq_ss_aln_j
		
		calc<-matrix(0, ncol=2, nrow= nb_aa, dimnames = list(names,c("r", "m")))
		calc[,2]<-trunc(mf) #inferior mf round
		calc[,1]<-(mf-calc[,2]) #remainder

		calc<-calc[order(calc[,1],rownames(calc), decreasing =T),] #tri des m par ordre decr des r puis pas ordre alpha inverse des aa
		
		#the sum of the "m" column must be equal to the sum of "n"
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
	
	if(double_passing){
	    # Calculating ELSC score for each valid position
	    for(i in nb_valid_pos:1){
		pos_i <- Valid_pos[i] #current valid position
		
		cat(paste("pos_i : ", pos_i, "\n"))
		
		for(j in nb_valid_pos:i){
		    pos_j <- Valid_pos[j]
		    
		    aln_tot_i<-AA[[pos_i]]
		    aln_tot_i[is.na(aln_tot_i)] <- 0
		    
		    aln_tot_j<-AA[[pos_j]]
		    aln_tot_j[is.na(aln_tot_j)] <- 0
		    
		    Nyj<-colSums(aln_tot_j) #number of amino acid y in position j in the full alignment
		    v<-sort(colSums(aln_tot_i),decreasing = TRUE)
		    v<-unique(c(which(aln_tot_i[,names(v[1])] == 1))) #sequences where appears the most occuring amino acid in position i

		    ss_aln_i<-aln_tot_i[v,]
		    nb_seq_ss_aln_i <- length(ss_aln_i[,1])
		    
		    ss_aln_j<-aln_tot_j[v,]
		    nb_seq_ss_aln_j <- length(ss_aln_j[,1])
		    
		    n<-colSums(ss_aln_j) #number of amino acid y in position j in the sub-alignment
		    mf<-(Nyj/nb_seq)*nb_seq_ss_aln_j
		    
		    calc<-matrix(0, ncol=2, nrow= nb_aa, dimnames = list(names,c("r", "m")))
		    calc[,2]<-trunc(mf) #inferior mf round
		    calc[,1]<-(mf-calc[,2]) #remainder

		    calc<-calc[order(calc[,1],rownames(calc), decreasing =T),] #tri des m par ordre decr des r puis pas ordre alpha inverse des aa
		    
		    #the sum of the "m" column must be equal to the sum of "n"
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
		    COV2[pos_i, pos_j]<- COV2[pos_i, pos_j] + (-log(prod(choose(Nyj, n)/choose(Nyj, m)))) #matrice asymetrique - vrais scores obtenus 
		}
	    }
	    
	    #Data is calculated twice then divided to get the mean
	    COV2 <- COV2/2
	}
	
	#Complete matrix second triangular part
	COV2 <- COV2 + t(COV2)
	diag(COV2) <- diag

	COV3<-(COV2+t(COV2))/2	
	COV4<-matrix(COV2,nrow=nb_pos)
	COV5<-matrix(COV2,nrow=nb_pos)
	for (i in 1:nb_pos) {
		for (j in i:nb_pos) {
			COV4[j,i]=COV2[i,j] 
			COV5[i,j]=COV2[j,i] 

		}
	}
	
	rownames(COV3)<-paste(rep("pos", nb_pos),colnames(N),sep="")
	colnames(COV3)<-paste(rep("pos", nb_pos),colnames(N),sep="")
	rownames(COV4)<-paste(rep("pos", nb_pos),colnames(N),sep="")
	colnames(COV4)<-paste(rep("pos", nb_pos),colnames(N),sep="")
	rownames(COV5)<-paste(rep("pos", nb_pos),colnames(N),sep="")
	colnames(COV5)<-paste(rep("pos", nb_pos),colnames(N),sep="")

	if(!is.null(fileCSV)){
		write.table(COV2, fileCSV, sep = "\t", quote = FALSE, col.names = NA)
	}

	#normalisation
# 	COV2<-1-(COV2/max(COV2))
# 	diag(COV2)<-0

# 	COV4<-1-(COV4/max(COV4))
# 	diag(COV4)<-0

	#Reducting the final correlation matrix to the valid positions
	COV2 <- COV2[Valid_pos,]
	COV2 <- COV2[,Valid_pos]
	
	#Reducting the final correlation matrix to the valid positions
	COV3 <- COV3[Valid_pos,]
	COV3 <- COV3[,Valid_pos]
	
	#Reducting the final correlation matrix to the valid positions
	COV4 <- COV4[Valid_pos,]
	COV4 <- COV4[,Valid_pos]

	#Removing "Inf" and "NaN" values
	COV2[is.infinite(COV2)] <- 0
	COV2[is.na(COV2)] <- 0
	
	#Removing "Inf" and "NaN" values
	COV3[is.infinite(COV3)] <- 0
	COV3[is.na(COV3)] <- 0
	
	#Removing "Inf" and "NaN" values
	COV4[is.infinite(COV4)] <- 0
	COV4[is.na(COV4)] <- 0
	
	
	res <- list()
	res$gross <- COV2
	
	if(z_score){
	    COV2 <- (COV2-mean(COV2))/sd(COV2)
	    COV3 <- (COV3-mean(COV3))/sd(COV3)
	    COV4 <- (COV4-mean(COV4))/sd(COV4)
	    COV5 <- (COV5-mean(COV5))/sd(COV5)
	}
	
	res$normalized <- COV2
	
	return(res)
}
