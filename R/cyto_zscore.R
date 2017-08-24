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
cyto_zscore <- function(corr_struct, contact= NULL, filepath){
  
  pos_names <- colnames(corr_struct)
  nb_pos <- length(pos_names)
  
  if(is.null(contact) | missing(contact)){
    positions <- pos_names
  
    write("pos_i (pp) pos_j = z_score", filepath, append= FALSE)
  
    for(pos_i in positions){
      j_positions <- positions[-1]
      positions <- positions[-1]
      
      for(pos_j in j_positions){
	val <- corr_struct[pos_i, pos_j]
	
	current_line <- paste(pos_i, "(pp)", pos_j, "=", val)
	write(current_line, filepath, append= TRUE)
      }
    }
  } else {
    i_positions <- contact$pairs_i
    j_positions <- contact$pairs_j
    nb_pairs <- length(i_positions)
    
    for(k in 1:nb_pairs){
      pos_i <- i_positions[k]
      pos_j <- j_positions[k]
      val <- corr_struct[pos_i, pos_j]
      
      current_line <- paste(pos_i, "(pp)", pos_j, "=", val)
      write(current_line, filepath, append= TRUE)
    }
    
  }
}
