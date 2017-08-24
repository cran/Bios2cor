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

entropy <- function (align, fileHelix= NULL) {
	nb_seq <- length(align)
	seq_size <- length(align[[1]])
	
	#base used in the log function (it refers to the number of amino acids)
	log_base <- 20

	#Apply functions version 
	entropy <- matrix(0, ncol= seq_size, nrow= 1)
	entropy <- unlist(lapply(1:seq_size, function(i){
	  checked_acids <- c()
	  entropy[i] <- -sum(unlist(lapply(1:nb_seq, function(j){
	    current_acid <- align[[j]][[i]]
	    if((current_acid != "-") & (!(current_acid %in% checked_acids))){
	      #occurrence of the current acid in the position j in all sequences
	      occ_current_acid <- sum(unlist(lapply(align, "[[", i)) == current_acid)
	      
	      #frequence of the acid in this position
	      freq_acid <- occ_current_acid/nb_seq
	      freq_acid_log <- freq_acid * log(freq_acid, base= log_base)
	      
	      #current acid has now been checked
	      checked_acids <<- c(checked_acids, current_acid)
	      
	      #adding the frequence of the current acid to the entropy score of this position
	      freq_acid_log
	    } else {
	      0
	    }
	  })))
	}))
	
	
	#Classic version (without apply functions)
# 	entropy <- matrix(0, ncol= seq_size, nrow= 1)
# 	for (i in 1:seq_size) {
# 		#list that contains the acids that have been counted while exploring the column
# 		checked_acids <- c()
# 		entropy[i] <- 0
# 		for (j in 1:nb_seq) {
# 			#name of the current acid
# 			current_acid <- align[[j]][[i]]
# 			
# 			#if current acid has not been checked yet and is different than "-"
# 			if((current_acid != "-") & (!(current_acid %in% checked_acids))){
# 				
# 				#occurrence of the current acid in the position j in all sequences
# 				occ_current_acid <- sum(unlist(lapply(align, "[[", i)) == current_acid)
# 				
# 				#frequence of the acid in this position
# 				freq_acid <- occ_current_acid/nb_seq
# 				freq_acid_log <- freq_acid * log(freq_acid, base= 20)
# 				
# 				#adding the frequence of the current acid to the entropy score of this position
# 				entropy[i] <- entropy[i] + freq_acid_log
# 				
# 				#current acid has now been checked
# 				checked_acids <- c(checked_acids, current_acid)
# 			}
# 		}
# 		#applying the minus one to each entropy value
# 		entropy[i] <- -entropy[i]
# 	}

	##Considering Ballesteros - Weinstein notation
	if(!is.null(fileHelix)){
	  #each line of the helix file looks like : "name1:x1-x2-x3"
	  Helix<-read.table(fileHelix, sep= "-", stringsAsFactors= FALSE)
	  nb_helix <- length(Helix[,1])
	  
	  #so we split each line like this : "name1:x1" "x2" "x3" (read.table can only have one separator)
	  #then we split the "name1:x1" : "name1" "x1" for each line of our data frame Helix we keep "x1"
	  lapply(1:nb_helix, function(x){ print(Helix[x,1]); Helix[x,1] <<- strsplit(Helix[x,1], split= ":")[[1]][2] })
	  
	  no_helix_inds <- unlist(lapply(1:(nb_helix-1), function(i){ if(as.numeric(Helix[i+1,1])-as.numeric(Helix[i,3]) > 1) c((as.numeric(Helix[i,3])+1):(as.numeric(Helix[i+1,1])-1)) }))
	  if(as.numeric(seq_size)-as.numeric(Helix[nb_helix,3]) > 1) no_helix_inds <- c(no_helix_inds, (as.numeric(Helix[nb_helix,3])+1):seq_size)
	  
	  helix_inds <- setdiff(1:seq_size, no_helix_inds)
	  
	  full_names <- unlist(lapply(1:7,function(i){paste(rep(i,length(Helix[i,1]:Helix[i,3])),(Helix[i,1]:Helix[i,3]+50-Helix[i,2]), sep=".")}))
	  
	  names <- c(1:seq_size)
	  
	  names[helix_inds] <- full_names
	} else {
	  names <- c(1:seq_size)
	}
	
	names <- paste("pos", names, sep="")
	names(entropy) <- names
	
	#Removing 0 entropy values
# 	entropy <- entropy[entropy != 0]
	
	return(entropy)
}
