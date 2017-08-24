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
entropy_graph <- function(entropy, corr, filter= NULL, elite, high, csv= TRUE, name){
  noms <- colnames(corr)
  taille <- length(noms)
  if(!is.null(filter)){
    for(i in 1:taille){
      for(j in 1:taille){
	nom_i <- noms[i]
	nom_j <- noms[j]
	ponderation_i <- filter[nom_i]
	ponderation_j <- filter[nom_j]
	corr[i,j] <- corr[i,j]*ponderation_i*ponderation_j
      }
    }
  }
  
  #Checking if ballesteros notation is used
  ballesteros <- grep(".", names(entropy)[1], fixed= TRUE)
  
  #Total number os positions
  nb_pos <- length(entropy)
  
  #Calculating possible entropy combinations
  i_inds <- unlist(lapply(1:nb_pos, function(x){ rep(x, nb_pos-x) }))
  j_inds <- unlist(lapply(2:nb_pos, function(x){ x:nb_pos }))
  
  #Entropy values
  i <- entropy[i_inds]
  i_names <- names(entropy)[i_inds]
#   print("i_names : ") ; print(i_names)
  
  j <- entropy[j_inds]
  j_names <- names(entropy)[j_inds]
#   print("j_names : ") ; print(j_names)
  
  #Combination numbers
  nb_comb <- length(i)
  
# A five-level color code indicativ of its rank (
# blue: the top 25 Z-scores,
# light blue: the top 275 Z-scores,
# pink: the bottom 275 Z-scores,
# red: the bottom 25 Z-scores,
# others: grey)
  
  #Collecting correlation score for each positions combination (excluding)
  if(length(ballesteros) != 0){
    
    #Positions with gaps could have been removed from correlation structure (corr)
#     unvalid_position_names <- setdiff(names(entropy), colnames(corr))
    
    corr_score <- c()
    for(k in 1:nb_comb){
      i_name <- i_names[k]
      j_name <- j_names[k]
      
      if(i_name %in% noms && j_name %in% noms){
	corr_score[k] <- corr[i_name, j_name]
      } else {
	corr_score[k] <- 0
      }
    }
  } else {
    print("no ballesteros detected")
    
    pos_names <- colnames(corr)
    nb_corr_pos <- length(pos_names)
    pos_helix <- unlist(lapply(1:nb_corr_pos, function(x){substr(pos_names[x], 4, nchar(pos_names[x]))})) #removing "pos" in each "pos36"...etc
    
    corr_score <- c()
    for(k in 1:nb_comb){
      i_name <- i_names[k]
      j_name <- j_names[k]
      
      if(i_name %in% noms && j_name %in% noms){
	corr_score[k] <- corr[i_name, j_name]
      } else {
	corr_score[k] <- 0
      }
    }
  }
  
  #Sorting considered positions
  sorted_corr_score <- sort(corr_score[corr_score != 0])
#   print("head(sorted_corr_score) = ") ; print(sorted_corr_score[1:40])
#   print("tail(sorted_corr_score) = ") ; print(tail(sorted_corr_score))
  
  nb_sorted_z_score <- length(sorted_corr_score)
  bound_1 <- sorted_corr_score[1]
  bound_2 <- sorted_corr_score[elite]
  bound_3 <- sorted_corr_score[nb_sorted_z_score-high]
  bound_4 <- sorted_corr_score[nb_sorted_z_score-elite]
  bound_5 <- sorted_corr_score[high]
  last_bound <- sorted_corr_score[nb_sorted_z_score]
  
#   print(paste("borne 1 : ", bound_1, sep= ""))
#   print(paste("borne 2 : ", bound_2, sep= ""))
#   print(paste("borne 3 : ", bound_3, sep= ""))
#   print(paste("borne 4 : ", bound_4, sep= ""))
#   print(paste("borne 5 : ", bound_5, sep= ""))
#   print(paste("borne 6 : ", last_bound, sep= ""))
  
  #Assigning points colors
  color_labels <- c("red", "pink", "lightblue", "blue", adjustcolor("grey", alpha.f= 0.3))
  neutral_color <- color_labels[5]
  my_colors <- c(rep(neutral_color, nb_comb))
  my_colors[corr_score >= bound_2 & corr_score <= bound_5] <- color_labels[2]
  my_colors[corr_score >= bound_3 & corr_score <= bound_4] <- color_labels[3]
  my_colors[corr_score >= bound_1 & corr_score <= bound_2] <- color_labels[1]
  my_colors[corr_score >= bound_4 & corr_score <= last_bound] <- color_labels[4]
  
  #Adapting points size
  points_size <- c(rep(1, nb_comb))
  points_size[my_colors == neutral_color] <- 1/3
  
  pdf(file= name)
    plot(i, j, main= basename(name), xlab= "S[i]", ylab= "S[j]", pch= 20, col= my_colors, cex= points_size)
  dev.off()
  
  if(csv){
    z_score <- scale(corr_score)
    name_csv <- paste(strsplit(name, "\\.")[[1]][1], ".csv", sep= "")
    write.table(list(i_names, j_names, corr_score, z_score, i, j), file= name_csv, row.names= FALSE, col.names= c("i_names", "j_names", "corr_score", "z_score", "entropy_i", "entropy_j"))
  }
}

