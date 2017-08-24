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
create_pcafile <- function(pca_struct, entropy= NULL, filepath, top_positions= 10, top_pca1= NULL, top_pca2= NULL){
  pca_coord_initial <- pca_struct$coord
  
  pca_positions <- rownames(pca_coord_initial)
  
  #Ordering data with first component
  ordered_res_1 <- order(pca_coord_initial[,1], decreasing= TRUE)
  pca_coord <- pca_coord_initial[ordered_res_1,]
  
  first_val <- abs(pca_coord[1,1])
  last_val <- abs(pca_coord[length(pca_positions), 1])
  if(first_val < last_val){
    ordered_res_1 <- order(pca_coord_initial[,1], decreasing= FALSE)
    pca_coord <- pca_coord_initial[ordered_res_1,]
  }
  
  pca_positions <- pca_positions[ordered_res_1]
  pca_size <- length(pca_coord[,1])
  pca_dim <- length(pca_coord[1,])
  
  head <- "position"
  if(!missing(entropy)){
    head <- paste(head, "entropy")
  }
  
  lapply(1:pca_dim, function(dim){
    head <<- paste(head, paste("PCA", dim, sep= ""))
  })
  
  write(head, file= filepath, append= FALSE)
  
  if(!missing(entropy)){
    for(pos in 1:pca_size){
      pos_line <- pca_coord[pos,]
      position <- pca_positions[pos]
      entropy_val <- entropy[position]
      
      #Ignoring possible NaN values
      if(sum(is.na(pos_line)) <= 0){
	coord_tmp <- paste(pca_coord[pos, 1:pca_dim], collapse= " ")
	current_line <- paste(position, entropy_val, coord_tmp)
	write(current_line, file= filepath, append= TRUE)
      }
    }
  } else {
    for(pos in 1:pca_size){
      pos_line <- pca_coord[pos,]
      position <- pca_positions[pos]
      
      #Ignoring possible NaN values
      if(sum(is.na(pos_line)) <= 0){
	coords <- ""
	for(dim in 1:pca_dim){
	  coords <- paste(coords, pca_coord[pos,dim])
	}
	current_line <- paste(position, coords, sep= "")
	write(current_line, file= filepath, append= TRUE)
      }
    }
  }
  
  ##Storing best pca 1 values file
  wanted_pca <- 2
  if(!missing(top_pca1)){
    if(!is.null(top_pca1)){
      if(pca_size < top_positions) top_positions <- pca_size
      header_pca_1 <- paste(c("position", paste("PCA_", 1:wanted_pca, sep= "")), collapse= " ")
      write(header_pca_1, file= top_pca1, append= FALSE)
      
      for(i in 1:top_positions){
	position <- pca_positions[i]
	coords <- paste(pca_coord[i, 1:wanted_pca])
	current_line <- paste(c(position, coords), collapse= " ")
	
	write(current_line, file= top_pca1, append= TRUE)
      }
    }
  }
  
  ##Storing best pca 2 values file
  #Ordering data with second component
  ordered_res_2 <- order(pca_coord_initial[,2], decreasing= TRUE)
  pca_coord <- pca_coord_initial[ordered_res_2,]
  pca_positions <- pca_positions[ordered_res_2]
  
  first_val <- abs(pca_coord[1,2])
  last_val <- abs(pca_coord[length(pca_positions), 2])
  if(first_val < last_val){
    ordered_res_2 <- order(pca_coord_initial[,2], decreasing= FALSE)
    pca_coord <- pca_coord_initial[ordered_res_2,]
  }
  
  pca_positions <- pca_positions[ordered_res_2]
  pca_size <- length(pca_coord[,1])
  pca_dim <- length(pca_coord[1,])
  
  wanted_pca <- 2
  if(!missing(top_pca2)){
    if(!is.null(top_pca2)){
      if(pca_size < top_positions) top_positions <- pca_size
      header_pca_2 <- paste(c("position", paste("PCA_", 1:wanted_pca, sep= "")), collapse= " ")
      write(header_pca_2, file= top_pca2, append= FALSE)
      
      for(i in 1:top_positions){
	position <- pca_positions[i]
	coords <- paste(pca_coord[i, 1:wanted_pca])
	current_line <- paste(c(position, coords), collapse= " ")
	
	write(current_line, file= top_pca2, append= TRUE)
      }
    }
  }
}

 
