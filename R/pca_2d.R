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
pca_2d <- function(pca_struct, abs= 1, ord= 2, filepath){
  if(is.null(abs)) abs <- 1
  if(is.null(ord)) ord <- 2
  
  pca_coord <- pca_struct$coord
  pca_abs <- pca_coord[, abs]
  pca_ord <- pca_coord[, ord]
  
  pca_x <- paste("PCA", abs, sep= "")
  pca_y <- paste("PCA", ord, sep= "")
  
  png(filepath)
    plot(pca_abs, pca_ord, xlab= pca_x, ylab= pca_y, main= basename(filepath))
  dev.off()
}
