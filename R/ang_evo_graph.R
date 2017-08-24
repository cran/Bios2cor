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
ang_evo_graph <- function(dynamic_struct, top_positions, filepath){
  #Importing rotamer informations
  tor <- dynamic_struct$tor
  
  nb_max_figures <- 35
  nb_col <- 5
  
  ##Generating ordered plots
  final_filepath <- filepath
  
  nb_top_positions <- length(top_positions)
  if(nb_top_positions > 20){
    fp_no_ext <- substr(filepath, 1, nchar(filepath)-3)
    final_filepath <- paste(fp_no_ext, "%03d", ".pdf", sep= "")
  }
  pdf(final_filepath)
    layout(matrix(1:nb_max_figures, ncol= nb_col, byrow= TRUE))
    
    #Setting special parameters
    par(mai= c(0.25, 0.25, 0.2, 0.2))
    
    #Alphanumerical graph ordering 
    graph_order <- c(1:nb_top_positions)
    graph_order <- sort(unlist(as.integer(lapply(strsplit(top_positions, "\\."), "[[", 1))), index.return= TRUE)$ix
    
    for(i in graph_order){
      position_name <- top_positions[i]
      position_angles <- tor[, position_name]
      nb_angles <- length(position_angles)
      plot(x= 1:nb_angles, y= position_angles, cex= 0.1, xlab= "frames", ylab="dihedral angles", main= position_name)
    }
    
  dev.off()
}
