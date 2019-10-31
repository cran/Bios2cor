#  
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
write.entropy <- function(entropy, filepathroot=NULL){

    if (missing(entropy)) {
      stop("A matrix of type 'entropy' is required.")  
    }

    if (is.null(filepathroot)) {
       filename <- "ENTROPY.csv"
    } else {
       filename <- paste(filepathroot,"_ENTROPY.csv",sep = "")
    }

    nb_pos <- length(entropy)
    positions <- names(entropy)
 
    head <- paste("position", "entropy")
    write(head, file= filename, append= FALSE)
    
    for(i in 1:nb_pos){
      pos <- positions[i]
	if (entropy[i]==''|is.na(entropy[i])){
      		val <- ''
        	current_line <- paste(pos, val)
        	write(current_line, file= filename, append= TRUE)
         } else {
     		val <- as.numeric(entropy[i])
                val <- format(val, digits=3, nsmall=3)
		current_line <- paste(pos, val)
      		write(current_line, file= filename, append= TRUE)
    	}
   }

  entropy <- as.numeric(na.omit(entropy))
    if (is.null(filepathroot)) {
       filename <- "ENTROPY_HIST.png"
    } else {
       filename <- paste(filepathroot,"_ENTROPY_HIST.png",sep = "")
    }

  png(filename)
    hist(entropy)
  dev.off()


}

 
