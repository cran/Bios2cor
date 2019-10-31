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
top_pairs_analysis <- function(corr_matrix, filepathroot= NULL, top= 25, entropy = NULL, filter = NULL ){
  
  if(missing(corr_matrix)){
    stop("A correlation matrix is required")
  }

  if(!is.matrix(corr_matrix)){
     stop("The first argument must be a score or Zscore matrix.")
  }

    if(is.null(filepathroot)){
      filename_contact <- paste("TOP", top, "_CONTACTS.csv", sep="")
      filename_score <- paste("TOP", top, "_SCORES.csv", sep="")
    } else {
      filename_contact <- paste(filepathroot,"_TOP", top, "_CONTACTS.csv", sep="")
      filename_score <- paste(filepathroot,"_TOP", top, "_SCORES.csv", sep="")
    }

 names <- colnames(corr_matrix)
  corr_matrix[lower.tri(corr_matrix)] <- 0
  

  corr_score <- corr_matrix
  names <- colnames(corr_score)
  size_matrix <- length(names)

  # Taking into account the delta filter 
  if(!is.null(filter)){
    for(i in 1:size_matrix){
      for(j in 1:size_matrix){
	name_i <- names[i]
	name_j <- names[j]
	ponderation_i <- filter[name_i]
	ponderation_j <- filter[name_j]
	# score equal to zero when at least one of the two positions is out of the allowed entropy range
	corr_score[i,j] <- corr_score[i,j]*ponderation_i*ponderation_j
      }
    }
  }



  # Extracting top value positions
  x <- which(corr_score >= sort(corr_score, decreasing= TRUE)[top], arr.ind= TRUE)
  
  nb_pairs <- length(x[,1])
  positions <- unique(c(x[,1], x[,2]))
  
  nb_positions <- length(positions)
  
  contacts <- c()
  for(i in 1:nb_positions){

    pos <- positions[i]
    nb_contacts <- 0
    for(j in 1:nb_pairs){
      node_1 <- x[j,1]
      node_2 <- x[j,2]
      if((pos == node_1) | (pos == node_2)){
	nb_contacts <- nb_contacts+1
      }
    }
    contacts[i] <- nb_contacts
  }
  
  # Writing results
  positions <- names[positions]
  
   if (is.null(entropy)){
 # Writing list of top positions with number of contacts 
    head <- paste("position", "contact")
    write(head, filename_contact, append= FALSE)
    for(k in 1:nb_positions){
      position <- positions[k]
      contact <- contacts[k]
      current_line <- paste(position, contact)
      write(current_line, filename_contact, append= TRUE)
    }
 } else {
   entropy <- entropy
   head <- paste("position", "contact", "entropy")
    write(head, filename_contact, append= FALSE)
    for(k in 1:nb_positions){
      position <- positions[k]
      contact <- contacts[k]
      entropy_k <- format(as.numeric(entropy[position]), digits=3, nsmall=3)
      current_line <- paste(position, contact, entropy_k)
      write(current_line, filename_contact, append= TRUE)
    }
}

  
  pairs_i <- x[,1] #pair_i
  pairs_j <- x[,2] #pair_j
  
  res <- list()
  res$pairs_i <- names[pairs_i]
  res$pairs_j <- names[pairs_j]
  res$positions <- positions
  res$contacts <- contacts
  
  ## Writing top score positions
    score <- unlist(lapply(1:nb_pairs, function(a){corr_score[x[a,1], x[a,2]]}))
    index <- sort(score, decreasing= TRUE, index.return= TRUE)$ix
    ordered_x <- x[index,]
    
    write("pair_i pair_j score_ij", filename_score, append= FALSE)
    for(i in 1:nb_pairs){
      pair_i <- names[ordered_x[i,1]]
      pair_j <- names[ordered_x[i,2]]
      score_ij <- score[index[i]]
      score_ij <- format(score_ij,digits=3, nsmall=3)
      score_line <- paste(pair_i, pair_j, score_ij)
      write(score_line, filename_score, append= TRUE)
    }


  
  return (res)
}
 
