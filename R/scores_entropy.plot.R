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
scores_entropy.plot <- function(entropy, corr_matrix, filepathroot=NULL, elite=25, high=275, filter=NULL){


  if(missing(entropy)) 
    stop("An 'entropy' object is required") 
  

  if(missing(corr_matrix)) 
    stop("A correlation matrix is required") 

  if(!is.matrix(corr_matrix)) 
    stop("corr_matrix must be a matrix. Select a score or Zscore matrix") 

  corr_score <- corr_matrix
  names_c <- colnames(corr_matrix)
  size_matrix <- length(names_c)
  names_s  <- names(entropy)


  if(!is.null(filter)){
    for(i in 1:size_matrix){
      for(j in 1:size_matrix){
	name_i <- names_c[i]
	name_j <- names_c[j]
	ponderation_i <- filter[name_i]
	ponderation_j <- filter[name_j]
	corr_score[i,j] <- corr_score[i,j]*ponderation_i*ponderation_j
      }
    }
  }
  
# Checking the numbering used in the entropy object and in the correlation matrix 
  numbering1 <- grep(".", names_s, fixed= TRUE)
  numbering2 <- grep(".", names_c, fixed= TRUE)
	
# Analysis is possible only if the same numbering is used for entropy and correlation
  if(((length(numbering1) == 0) && (length(numbering2) == 0)) || ((length(numbering1) != 0) && (length(numbering2) != 0))){

     entropy_i <- c()
     entropy_j <- c()
     score <- c()

      k <-0
      for(i in 1:(size_matrix-1)){
         for(j in (i+1):size_matrix){
	  posi <- names_c[i]
	  posj <- names_c[j]
          #find entropy for each position of the pair
	  if(posi %in% names_s && posj %in% names_s){
	     k <- k+1
	     entropy_i[k] <- entropy[posi]
	     entropy_j[k] <- entropy[posj]
	     score[k] <- corr_score[posi, posj]
	  }else{
	     next
	  }
        }
      }


results <- cbind(score,entropy_i,entropy_j)
result_ordered <- results[order(as.numeric(results[,1]), decreasing = TRUE),]

x<- result_ordered[,1]
i<-result_ordered[,2]
j<-result_ordered[,3]

nb_pairs <-length(x)

bottom_elite <- nb_pairs - elite
bottom_high <- nb_pairs - high
stop_high <- bottom_elite - 1
start_high <- elite + 1

if (is.null(filepathroot)){
   filename <- "ei_ej.pdf"
} else {
filename = paste(filepathroot, "_ei_ej.pdf", sep ="")
}

pdf(file= filename)

plot(i, j, main= filepathroot, xlab= "S[i]", ylab= "S[j]", pch= ".", col= "gray80",cex.axis=1.4, cex.lab=1.4)
if(is.null(filter)) {
points(i[bottom_high:stop_high], j[bottom_high:stop_high], col = "lightpink1", pch=20)
points(i[bottom_elite:nb_pairs], j[bottom_elite:nb_pairs], col = "red", pch=20)
points(i[start_high:high], j[start_high:high], col = "skyblue2", pch=20)
points(i[1:elite], j[1:elite], col = "blue", pch=20)
} else {
points(i[start_high:high], j[start_high:high], col = "skyblue2", pch=20)
points(i[1:elite], j[1:elite], col = "blue", pch=20)
}

dev.off()

  } else {
 print("Mismatch in the notation used for the correlation and entropy files. Please verify your files!") 
 }
}
