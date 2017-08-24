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
centered_pca <- function(mat, m = NULL, pc= 3, dec_val= 3, eigenvalues_csv= NULL){
  
  #Centering the initial matrix
  cor <- mat
  size <- length(mat[,1])
  
  #Positions count
  if(is.null(pc) | missing(pc)) pc <- size
  print(paste("pc : ", size))
  
  #identity matrix
  I <- diag(1, size)

  #mat matrix of ones
  ONES <- matrix(1, nrow = size, ncol = 1)
  
  #Without mentioning any filter, every position has the same weighting
  if(is.null(m) | missing(m)){
    print("No filter applied")
    m <- matrix(1/size, nrow= size, ncol= 1)
    mat_names <- colnames(mat) #position names of mat
  } else {
    print("Filter applied")
    #Adapting the m filter to contain the positions of mat (positions that have a correct gap proportion)
    mat_names <- colnames(mat) #position names of mat
    
    m <- m[mat_names]
  }
  
  #m vector of mass
  BigI<-I-(ONES%*%t(m))
  
  #compute mat cross-product matrix
  S <-  BigI %*% cor %*% t(BigI)
  
  #Diagonalizing the centered matrix
  res <- list()
  eigen <- eigen(S)
  
  res$eigen <- round(eigen$values[1:pc], 3)
  eigen.perc <- (abs(eigen$values) * 100) / sum(eigen$values[eigen$values>0])
  res$eigen.perc <- round(eigen.perc[1:pc], 3)
  
  all_eigen_values <- eigen$values
  nb_eigenvalues <- length(all_eigen_values)
  
  #only positive eigenvalues are kept
  eigen$vectors <- eigen$vectors[, eigen$values > 0]
  eigen$values <- eigen$values[eigen$values > 0]
  nb_positiv_eigenvalues <- length(eigen$values)
  
  #Storing eigen values in csv file
  if(!is.null(eigenvalues_csv)){
    positiv_ev <- sum(eigen$values)
    perc_positiv_ev <- unlist(lapply(eigen$values, function(x){x*100/positiv_ev}))
    
    csv_tab <- data.frame("eigen$values"= all_eigen_values, "positiv_percent"= c(perc_positiv_ev, rep(0, nb_eigenvalues-nb_positiv_eigenvalues)))
    write.table(csv_tab, row.names= FALSE, file= basename(eigenvalues_csv))
    
    eigenvalues_png <- paste(substr(eigenvalues_csv, 1, nchar(eigenvalues_csv)-4), ".png", sep= "")
    png(eigenvalues_png)
      plot(1:nb_eigenvalues, all_eigen_values, main= basename(eigenvalues_png))
      abline(h= 0, lty= 2)
    dev.off()
  }
  
  res$source<-list()
  res$source$cor <- cor
  res$source$m<-m
  
  #check principal components
  if (pc < 2) pc <- 3
  if (pc > length(eigen$values)) pc <- length(eigen$values)

  #Computing PCA for the centered and diagonalized matrix
  #eigenvalues are transformed into percentage
  #compute mat matrix of factor scores
  
  F <- diag(as.vector(m)^(-0.5)) %*% eigen$vectors %*% diag(eigen$values^0.5)
  
  coord <- data.frame(F[, 1:pc])
  rownames(coord) <- rownames(mat)
  colnames(coord) <- paste ("PC", (1:pc), sep = "")
  res$coord = round(coord, dec_val)
  class (res) <- c("pca")
  return (res)
}
