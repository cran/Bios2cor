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
centered_pca <- function(corr_matrix, filepathroot = NULL, filter = NULL, pc= NULL, dec_val= 5){
 
  if(missing(corr_matrix) | is.null(corr_matrix)) {
    stop("A correlation matrix is required")
  }

  if(!is.matrix(corr_matrix)){
     stop("The first argument must be a score or Zscore matrix.")
  }

 
 #Centering the initial matrix
  cor <- corr_matrix
  size <- length(corr_matrix[,1])
  
  #Positions count
  if(is.null(pc)) {
  pc <- size
  print(paste("pc : ", size))
  }

  #identity matrix
  I <- diag(1, size)

  #matrix of ones
  ONES <- matrix(1, nrow = size, ncol = 1)
  
  
  if(is.null(filter)) {
    print("No filter applied")
    #Without filter, the elements have the same mass/weight
    m <- matrix(1/size, nrow= size, ncol= 1)
    mat_names <- colnames(corr_matrix) # names of elements in corr_matrix
  } else {
    print("Filter applied")
    mat_names <- colnames(corr_matrix) #names of elements in corr_matrix
    #filter is limited to the elements present in corr_matrix (may differ in sequence analysis)
    m <- filter[mat_names]    
    SUM <- sum(m)
    m <- m/SUM
  }
  
  #m vector of mass
  BigI<-I-(ONES%*%t(m))
  
  #compute cross-product matrix
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
  
  #Storing and plotting eigen values 
    if(is.null(filepathroot)){
      eigen_csv <- "EIGEN.csv"
      eigen_png <- "EIGEN.png"
    } else {
      eigen_csv <- paste(filepathroot, "_EIGEN.csv", sep="")
      eigen_png <- paste(filepathroot, "_EIGEN.png", sep="")
    }
 
    positiv_ev <- sum(eigen$values)
    perc_positiv_ev <- unlist(lapply(eigen$values, function(x){x*100/positiv_ev}))
    
    png(eigen_png)
      plot(1:nb_eigenvalues, all_eigen_values, main= basename(eigen_png))
      abline(h= 0, lty= 2)
    dev.off()

    all_eigen_values <- round(all_eigen_values,digits = 3)
    perc_positiv_ev <- round(perc_positiv_ev, digits= 3)

    csv_tab <- data.frame("eigen$values"= all_eigen_values, "positiv_percent"= c(perc_positiv_ev, rep(0, nb_eigenvalues-nb_positiv_eigenvalues)))
    write.table(csv_tab, row.names= FALSE, file= eigen_csv)
    
 

  res$source<-list()
  res$source$cor <- cor
  res$source$m<-m
  
  #check principal components
  if (pc < 2) pc <- 3
  if (pc > length(eigen$values)) pc <- length(eigen$values)

  #compute the matrix of factor scores
  
  F <- diag(as.vector(m)^(-0.5)) %*% eigen$vectors %*% diag(eigen$values^0.5)
  
  coord <- data.frame(F[, 1:pc])
  rownames(coord) <- rownames(corr_matrix)
  colnames(coord) <- paste ("PC", (1:pc), sep = "")

  res$coord = round(coord, dec_val)
  class (res) <- c("pca")
  return (res)


 
}
