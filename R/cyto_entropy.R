# # Bioscor is free software: you can redistribute it and/or modify
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
cyto_entropy <- function(entropy, contact= NULL, filepath){
  
  pos_names <- names(entropy)
  nb_pos <- length(entropy)
  
  write("pos = entropy", filepath, append= FALSE)
  
  if(is.null(contact) | missing(contact)){
    positions <- pos_names
  } else {
    positions <- contact$positions
  }
  
  for(pos in positions){
    val <- entropy[pos]
    current_line <- paste(pos, "=", val)
    write(current_line, filepath, append= TRUE)
  }
}
