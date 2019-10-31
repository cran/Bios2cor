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
write.pca <- function(corr_pca, filepathroot=NULL, pc= NULL, entropy= NULL){
  
  if (missing(corr_pca)) {
      stop("A PCA object created by the centered_pca function is required")
  }

 if (is.null(filepathroot)) {
  filename <- "PCA_COORD.csv"
  }else{
  filename <-paste(filepathroot, "_PCA_COORD.csv", sep="")
  }


  pca_coord <- corr_pca$coord
  pca_positions <- rownames(pca_coord)
  pca_size <- length(pca_coord[,1])
  
  if (is.null(pc)) { 
  pca_dim <- length(pca_coord[1,])
  } else {
  pca_dim <- pc
  }

  head <- "position"
  if(!is.null(entropy)) {
    head <- paste(head, "entropy")
  }
  
  lapply(1:pca_dim, function(dim){
    head <<- paste(head, paste("PCA", dim, sep= ""))
  })
  
  write(head, file= filename, append= FALSE)
  
  if(!is.null(entropy)) {
    for(pos in 1:pca_size){
      pos_line <- pca_coord[pos,]
      position <- pca_positions[pos]
      entropy_val <- format(as.numeric(entropy[position]), digits=3, nsmall=3)
      
      #Ignoring possible NaN values
      if(sum(is.na(pos_line)) <= 0){
	coord_tmp <- paste(pca_coord[pos, 1:pca_dim], collapse= " ")
	current_line <- paste(position, entropy_val, coord_tmp)
	write(current_line, file= filename, append= TRUE)
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
	write(current_line, file= filename, append= TRUE)
      }
    }
  }
  
}

 
