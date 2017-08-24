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
shuffle_all <- function(align, method, fileHelix= NULL, gap_val= 0.8, z_score= TRUE, nb_iterations= 5){
nb_seq <- length(align)
nb_pos <- length(align[[1]])

#Initialization
corr <- NULL
entropy <- NULL

for(i in 1:nb_iterations){
  print(paste(method, "shuffle all : ", i, "/", nb_iterations))
  
  ##Correlation matrix part
  #randomizing align
  align_int <- random.msa(nb.seq= nb_seq, nb.pos= nb_pos, align= align)
  
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
  
  ##Entropy part
  entropy_int <- entropy(align, fileHelix= fileHelix)
  
  #adding intermediary entropy score to result
  if(is.null(entropy)){
    entropy <- entropy_int
  } else {
    entropy <- entropy + entropy_int
  }
}

#calculating correlation matrix mean
corr$gross <- corr$gross / nb_iterations
   if(z_score){
      corr$normalized <- corr$normalized / nb_iterations
   }
corr_position_names <- colnames(corr$gross)

#calculating entropy scores mean
entropy <- entropy[corr_position_names]
entropy <- entropy / nb_iterations

res <- list()
res$corr$gross<- corr$gross
if(z_score){
      res$corr$normalized <- corr$normalized
   }

res$entropy <- entropy

return (res)
}
