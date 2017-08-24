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
corr_contact <- function(corr, top= 25, contact_filepath= NULL, score_top_filepath= NULL){
  names <- colnames(corr)
  
  corr[lower.tri(corr)] <- 0
  
  ##Extracting top value positions
  x <- which(corr >= sort(corr, decreasing= TRUE)[top], arr.ind= TRUE)
  
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
  
  positions <- names[positions]
  
  if(!is.null(contact_filepath) & !missing(contact_filepath)){
    head <- paste("position", "contact")
    write(head, contact_filepath, append= FALSE)
    for(k in 1:nb_positions){
      position <- positions[k]
      contact <- contacts[k]
      current_line <- paste(position, "=", contact)
      write(current_line, contact_filepath, append= TRUE)
    }
  }
  
  pairs_i <- x[,1] #pair_i
  pairs_j <- x[,2] #pair_j
  
  res <- list()
  res$pairs_i <- names[pairs_i]
  res$pairs_j <- names[pairs_j]
  res$positions <- positions
  res$contacts <- contacts
  
  ##Score top positions
  if(!is.null(score_top_filepath) & !missing(score_top_filepath)){
    score <- unlist(lapply(1:nb_pairs, function(a){corr[x[a,1], x[a,2]]}))
    index <- sort(score, decreasing= TRUE, index.return= TRUE)$ix
    ordered_x <- x[index,]
    
    write("pair_i pair_j score_ij", score_top_filepath, append= FALSE)
    for(i in 1:nb_pairs){
      pair_i <- names[ordered_x[i,1]]
      pair_j <- names[ordered_x[i,2]]
      score_ij <- score[index[i]]
      
      score_line <- paste(pair_i, pair_j, score_ij)
      write(score_line, score_top_filepath, append= TRUE)
    }
  }
  
  return (res)
}
 
