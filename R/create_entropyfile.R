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
create_entropyfile <- function(entropy, filepath){
    nb_pos <- length(entropy)
    positions <- names(entropy)
    
    head <- paste("position", "entropy")
    write(head, file= filepath, append= FALSE)
    
    for(i in 1:nb_pos){
      pos <- positions[i]
      val <- entropy[i]
      
      #Creating the line to insert to the file
      current_line <- paste(pos, val)
      
      write(current_line, file= filepath, append= TRUE)
    }
}

 
