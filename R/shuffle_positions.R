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
shuffle_positions <- function(align, method, fileHelix= NULL, gap_val= 0.8, z_score= TRUE, nb_iterations= 5){
nb_seq <- length(align)
nb_pos <- length(align[[1]])

#Initialization
corr <- NULL

for(i in 1:nb_iterations){
  print(paste(method, "shuffle positions : ", i, "/", nb_iterations))
  
  #Transforming align element into matrix
  align_int <- sapply(1:nb_pos, function(x){unlist(lapply(align, "[[", x))})

  #Sampling align
  align_int <- sapply(1:nb_pos, function(x){align_int[,x] <- sample(align_int[,x])})
  
  #Rebuilding align structure
  align_int <- lapply(1:nb_seq, function(x){align_int[x,]})

  #calculating correlation matrix
  corr_int <- switch(method,
	"OMES"= omes(align_int, fileHelix= fileHelix, gap_val= gap_val, z_score= z_score),
	"MIP"= mip(align_int, fileHelix= fileHelix, gap_val= gap_val, z_score= z_score),
	"MCBASC"= mcbasc(align_int, fileHelix= fileHelix, gap_val= gap_val, z_score= z_score),
	"ELSC"= elsc(align_int, fileHelix= fileHelix, gap_val= gap_val, z_score= z_score)
  )
  
   #adding intermediary matrix to result
  if(is.null(corr)){
    corr$gross <- corr_int$gross
    if(z_score){
      corr$normalized <- corr_int$normalized
      }
  } else {
    corr$gross <- corr$gross + corr_int$gross
    if(z_score){
      corr$normalized <- corr$normalized + corr_int$normalized
      }
  }
  
}

#calculating correlation matrix mean
corr$gross <- corr$gross / nb_iterations
   if(z_score){
      corr$normalized <- corr$normalized / nb_iterations
   }

res <- list()

res$corr$gross<- corr$gross
if(z_score){
      res$corr$normalized <- corr$normalized
   }


return (res)
}
