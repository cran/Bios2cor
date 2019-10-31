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
scores.boxplot <- function(corr_matrix_list, name_list, filepathroot=NULL, elite=25, high=275){
  if (is.null(filepathroot)) {
    filename <- "BOXPLOT.png"
  }else{
    filename <- paste(filepathroot, "_BOXPLOT.png", sep = "")
  }
  corr_tab <- sapply(corr_matrix_list, function(x){ as.vector(x)})
  nb_objects <- length(corr_matrix_list)
  

  print(paste("Boxplot elements :", nb_objects))
  
  png(filename, width = 600, height = 400, units = "px", pointsize = 12)
    boxplot(corr_tab, names = name_list, xlab="Z-score", cex.lab=1.5, cex.axis=0.90, col = "grey", outcol = "grey", horizontal=T, las=1)
    for(i in 1:nb_objects){
      increasing_corr <- sort(corr_tab[,i])
      decreasing_corr <- sort(corr_tab[,i], decreasing= TRUE)
      
      top_high <- decreasing_corr[(elite+1):high]
      points(top_high, rep(i,length(top_high)), pch=16, col="dodgerblue")
      
      top_elite <- decreasing_corr[1:elite]
      points(top_elite, rep(i,length(top_elite)), pch=16, col="blue")
      
      bottom_high <- increasing_corr[(elite):high]
      points(bottom_high, rep(i,length(bottom_high)), pch=16, col="pink")
      
      bottom_elite <- increasing_corr[1:elite]
      points(bottom_elite, rep(i,length(bottom_elite)), pch=16, col="red")
    }
    
    #Draw a line at position v
    abline(v=4, lty=2)			
    
    #Box with special line width
    box(lwd = 3)
  dev.off()
}
