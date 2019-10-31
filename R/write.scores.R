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
write.scores <- function(correlation, entropy= NULL, filepathroot=NULL){
 
  if (missing(correlation)) {
     stop("A correlation object is required.")
  } 


  if (is.null(filepathroot)) {
     filename <- "CORR_SCORES.csv"
  }else{
     filename <-paste(filepathroot, "_CORR_SCORES.csv", sep="")
  }


  corr_score <- correlation$score
  corr_Zscore <- correlation$Zscore  

  names_c <- colnames(corr_score)
  nb_pos <- length(names_c)
  numbering1 <- grep(".", names_c, fixed= TRUE)
  
  if (!is.null(entropy)){
     names_s  <- names(entropy)
     numbering2 <- grep(".", names_s, fixed= TRUE)
     # Check that both entropy and correlation objects are from a MSA or a trajectory  
     if(((length(numbering1) == 0) && (length(numbering2) != 0)) || ((length(numbering1) != 0) && (length(numbering2) == 0))){
         print("Mismatch in the notation used for the correlation and entropy files. Please verify your files!") 
     }
  }	

 
   if(length(numbering1) != 0){                            
 
     # Read the score_noauto and the Zscore_noauto matrices
      corr_score_noauto <- correlation$score_noauto
      corr_Zscore_noauto <- correlation$Zscore_noauto  

      if (!is.null(entropy)){
         head <- paste("pos_i", "pos_j", "score", "Zscore", "score_noauto", "Zscore_noauto", "entropy_i", "entropy_j")
         write(head, file= filename, append= FALSE)
	  
	 for(i in 1:nb_pos){
            for(j in i:nb_pos){
	       pos_i <- names_c[i]
	       pos_j <- names_c[j]
	       if(pos_i != pos_j){
	  	  entropy_i <- format(as.numeric(entropy[pos_i]), digits=3, nsmall=3)
	  	  entropy_j <- format(as.numeric(entropy[pos_j]), digits=3, nsmall=3)
	  	  score <- format(as.numeric(corr_score[pos_i, pos_j]), digits=3, nsmall=3)
	  	  Zscore <- format(as.numeric(corr_Zscore[pos_i, pos_j]), digits=3, nsmall=3)
	          score_noauto <- format(as.numeric(corr_score_noauto[pos_i, pos_j]), digits=3, nsmall=3)
	          Zscore_noauto <- format(as.numeric(corr_Zscore_noauto[pos_i, pos_j]), digits=3, nsmall=3)
	  
	  	  # Create the line to insert to the file
	  	  current_line <- paste(pos_i, pos_j, score, Zscore, score_noauto, Zscore_noauto, entropy_i, entropy_j) 
	  	  write(current_line, file= filename, append= TRUE)
	       }
            }
         }
      } else {
         head <- paste("pos_i", "pos_j", "score", "Zscore", "score_noauto", "Zscore_noauto")
         write(head, file= filename, append= FALSE)
	  
	 for(i in 1:nb_pos){
            for(j in i:nb_pos){
	       pos_i <- names_c[i]
	       pos_j <- names_c[j]
	       if(pos_i != pos_j){
	  	  score <- format(as.numeric(corr_score[pos_i, pos_j]), digits=3, nsmall=3)
	  	  Zscore <- format(as.numeric(corr_Zscore[pos_i, pos_j]), digits=3, nsmall=3)
	          score_noauto <- format(as.numeric(corr_score_noauto[pos_i, pos_j]), digits=3, nsmall=3)
	          Zscore_noauto <- format(as.numeric(corr_Zscore_noauto[pos_i, pos_j]), digits=3, nsmall=3)
	  
	  	  # Create the line to insert to the file
	  	  current_line <- paste(pos_i, pos_j, score, Zscore, score_noauto, Zscore_noauto) 
	  	  write(current_line, file= filename, append= TRUE)
	       }
            }
         }
      }

   } else {

      if (!is.null(entropy)){
         head <- paste("pos_i", "pos_j", "score", "Zscore", "entropy_i", "entropy_j")
         write(head, file= filename, append= FALSE)

         for(i in 1:nb_pos){
            for(j in i:nb_pos){
	       pos_i <- names_c[i]
	       pos_j <- names_c[j]
	       if(pos_i != pos_j){
	  	  entropy_i <- entropy[pos_i]
	  	  entropy_j <- entropy[pos_j]
	  	  score <- format(as.numeric(corr_score[pos_i, pos_j]), digits=3, nsmall=3)
	   	  Zscore <- format(as.numeric(corr_Zscore[pos_i, pos_j]), digits=3, nsmall=3)
	          if (entropy_i == '') {
	 		entropy_i <- 'ND'
     		  }else{
	  	  entropy_i <- format(as.numeric(entropy_i), digits=3, nsmall=3)
		  }
	          if (entropy_j == '') {
	 		entropy_j <- 'ND'
     		  }else{
	  	  entropy_j <- format(as.numeric(entropy_j), digits=3, nsmall=3)
		  }
	  	  #Creating the line to insert to the file
	  	  current_line <- paste(pos_i, pos_j, score, Zscore, entropy_i, entropy_j)
	  	  write(current_line, file= filename, append= TRUE)
	       }
            }
         }
      }else{
         head <- paste("pos_i", "pos_j", "score", "Zscore")
         write(head, file= filename, append= FALSE)

   	 for(i in 1:nb_pos){
      	    for(j in i:nb_pos){
		pos_i <- names_c[i]
		pos_j <- names_c[j]
		if(pos_i != pos_j){
	  	   score <- format(as.numeric(corr_score[pos_i, pos_j]), digits=3,nsmall=3)
	  	   Zscore <- format(as.numeric(corr_Zscore[pos_i, pos_j]), digits=3,nsmall=3)
	  
	 	   #Creating the line to insert to the file
		   current_line <- paste(pos_i, pos_j, score, Zscore)
		   write(current_line, file= filename, append= TRUE)
	        }
      	     }
    	 }
      }
 
   } 

}

 
