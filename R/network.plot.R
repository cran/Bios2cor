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
network.plot <- function(top_pairs, filepathroot=NULL){
  pairs_i <- top_pairs$pairs_i
  pairs_j <- top_pairs$pairs_j

   if (is.null(filepathroot)) {
      filename <- "NETWORK.pdf"
   } else {
      filename <- paste(filepathroot, "_NETWORK.pdf", sep="")
   }

  pdf(filename)
    # create data:
    links <- data.frame(source= pairs_i, target= pairs_j)
    
    # Turn it into igraph object
    network <- graph_from_data_frame(d= links, directed= F) 
    
    # Count the number of degree for each node:
    deg <- degree(network, mode= "all")
    
    # Plot
    plot(network, vertex.size= deg*6, vertex.color= rgb(0.1,0.7,0.8,0.5))
  
  dev.off()
}
