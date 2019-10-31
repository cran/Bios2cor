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
pca_screeplot <- function(pca_struct, filepathroot=NULL){
  pca_coord <- pca_struct$coord
  
  nb_component <- length(pca_coord[1,])
  pca_positions <- rownames(pca_coord)

  components <- 1:nb_component
  
  variances <- unlist(lapply(1:nb_component, function(x){var(pca_coord[,x])}))

  if(is.null(filepathroot)) {
    filename <- "SCREEPLOT.png" 
  } else {
    filename <- paste(filepathroot, "_SCREEPLOT.png", sep="")
  }  
  
  png(filename, width = 600, height = 600, units = "px")
    plot(components, variances, main= basename(filename))
    lines(components, variances)
  dev.off()
}

 
