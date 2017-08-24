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
create_screeplot <- function(pca_struct, pca_index= NULL, filepath){
  pca_coord <- pca_struct$coord
  
  nb_component <- length(pca_coord[1,])
  pca_positions <- rownames(pca_coord)

  if(!is.null(pca_index)){
    pca_coord <- pca_coord[,pca_index]
    components <- pca_index
  } else {
    components <- 1:nb_component
  }
  
  variances <- unlist(lapply(1:nb_component, function(x){var(pca_coord[,x])}))
  
  png(filepath, width = 400, height = 800, units = "px")
    plot(components, variances, main= basename(filepath))
    lines(components, variances)
  dev.off()
}

 
