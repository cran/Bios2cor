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
angles.plot <- function(dynamic_structure, angles=NULL, filepathroot=NULL){
 
if (missing(dynamic_structure))
  stop("A 'dynamic_structure' object is required")

if (is.null(angles)){
  angles <- dynamic_structure$tor.names
  }


if (is.null(filepathroot)) {
  final_filepath <- "ANGLES.pdf"
  final_csv <- "ANGLES.csv" 
} else {
  final_filepath <- paste(filepathroot, "_ANGLES.pdf",sep="")
  final_csv <- paste(filepathroot, "_ANGLES.csv",sep="")
}


  #Importing torsional and angles data
  tor <- dynamic_structure$tor 
  nb_frames <- length(tor[,1])
  nb_angles <- length(angles)
  frames <- dynamic_structure$frames

  #Generating ordered plots
  nb_max_figures <- 35
  nb_col <- 5
  
    pdf(final_filepath)
    layout(matrix(1:nb_max_figures, ncol= nb_col, byrow= TRUE))     # plot up to 7 rows of 5 figures / page
    
    #Setting the margins of the graph in inches (mai)  
    par(mai= c(0.25, 0.25, 0.2, 0.2))
    
    #Alphanumerical graph ordering 
    graph_order <- c(1:nb_angles)
    graph_order <- sort(unlist(as.integer(lapply(strsplit(angles, "\\."), "[[", 1))), index.return= TRUE)$ix
    
    for(i in graph_order){
      angle_name <- angles[i]
      dihedral_angles <- tor[, angle_name]
      plot(x= 1:nb_frames, y= dihedral_angles, cex= 0.1, xlab= "frames", ylab="dihedral angles", main= angle_name)
    }
    
  dev.off()


 tor <- tor[,angles]    
 write.csv(cbind(frames,round(tor, digits=2)), file = final_csv)



}
