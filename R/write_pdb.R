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
write_pdb <- function(pca, trio_comp= c(1:3), pdb_filepath, pml_filepath){
  coords <- pca$coord
  pca_1 <- trio_comp[1]
  pca_2 <- trio_comp[2]
  pca_3 <- trio_comp[3]
  
  pos_names <- rownames(coords)
  nb_pos <- length(pos_names)
  nb_helix <- 7
  pos_prefix <- substr(pos_names[1], 1, 3) == "pos"
  
  #check if ballesteros notation is used or not
  ballesteros <- grep(".", pos_names, fixed= TRUE)
  nb_ballesteros <- length(ballesteros)
  
  any_ballesteros <- 0
  
  if(nb_ballesteros > 0) any_ballesteros <- 1
  
  if(any_ballesteros) {
    pos <- strsplit(pos_names, split= "[.]")
    
    #positions with true ballesteros notation
    ballesteros_positions <- unlist(lapply(pos, function(x){length(x) > 1}))
    
    #extracting positions helix numbers
    helix_numbers <- unlist(lapply(pos[ballesteros_positions], "[[", 1))
    if(pos_prefix) helix_numbers <- substr(helix_numbers, 4, nchar(helix_numbers))
    
    #positions without ballesteros notation have 0 as helix number
    pos_helix <- rep(0, nb_pos)
    pos_helix[ballesteros_positions] <- helix_numbers
    
    #extracting residue numbers
    res <- unlist(lapply(pos, "[[", 1))
    #removing "pos" in each element
    if(pos_prefix) res <- substr(res, 4, nchar(res))
    #pasting residue numbers of each position using ballesteros notation
    res[ballesteros_positions] <- unlist(lapply(pos[ballesteros_positions], "[[", 2))
    
    pos_residue <- res
    
    #letter conversion of helixes (in "1.36", "1" will become "A")
    alpha <- toupper(letters[1:nb_helix])
    
    pos_helix_letter <- rep("Z", nb_pos)
    pos_helix_letter[pos_helix != 0] <- alpha[as.integer(pos_helix[pos_helix != 0])]
    
#     print("pos_helix : ") ; print(pos_helix)
#     print("pos_helix_letter : ") ; print(pos_helix_letter)
#     print("pos_residue : ") ; print(pos_residue)
    
  } else {
    if(pos_prefix) pos_names <- unlist(lapply(1:nb_pos, function(x){substr(pos_names[x], 4, nchar(pos_names[x]))}))
    pos_residue <- pos_names
    pos_helix_letter <- rep("A", nb_pos)
  }
  
  if(any_ballesteros){
    colors <- c("blue", "red", "green", "yellow", "brown", "pink", "orange", "white")
    
    background_color <- paste("bg_color", colors[8])
    write(background_color, file= pml_filepath, append= TRUE)
  
    #Helix coloration
    for(i in 1:nb_helix){
      chain_coloration <- paste("color", colors[i], ", chain", alpha[i])
      write(chain_coloration, file= pml_filepath, append= TRUE)
    }
    
    #Configuring elements
    as_sphere <- "as sphere"
    sphere_scale <- "set sphere_scale, 0.3"
    write(as_sphere, file= pml_filepath, append= TRUE)
    write(sphere_scale, file= pml_filepath, append= TRUE)
  }
  
  ##fields required in each line
  #fields that doesn't change
  head <- "ATOM"
  atom_name <- "O"
  residue_name <- "HOH"
  occupation <- "1.00"
  B_factor <- 100.00
  chain_name <- "O"
  loc_indicator <- ""
  seg_ID <- ""
  charge <- ""
  iCode <- ""
  format <- "%-6s%5d %-4s%-1s%3s %-1s%4s%1s   %8s%8s%8s%6s%6s      %-4s%2s%2s"
  
  #fields that does change
  for(pos in 1:nb_pos){
    atom_number <- pos_names[pos]
    x <- coords[pos, pca_1]
    y <- coords[pos, pca_2]
    z <- coords[pos, pca_3]
    chain_ID <- pos_helix_letter[pos]
    residue_number <- pos_residue[pos] #ballesteros notation
    
    pdb_line <- sprintf(format, head, pos, atom_name, loc_indicator, residue_name, chain_ID, residue_number, iCode, x, y, z, occupation, B_factor, seg_ID, chain_name, charge)
    
    write(pdb_line, file= pdb_filepath, append= TRUE)
  }
  
  #Creating elements used to create axis
  max_x <- max(abs(coords[,1]))
  norm <- max_x*2
  axis_fake_atoms <- c("O","C","N", "Po")
  axis_fake_chain_ID <- "Z"
  x_axis_ID <- nb_pos+1
  y_axis_ID <- nb_pos+2
  z_axis_ID <- nb_pos+3
  origin_ID <- nb_pos+4
  
  x_axis <- sprintf(format, head, x_axis_ID, axis_fake_atoms[1], loc_indicator, residue_name, axis_fake_chain_ID, x_axis_ID, iCode, norm, 0, 0, occupation, B_factor, seg_ID, axis_fake_atoms[1], charge)
  y_axis <- sprintf(format, head, y_axis_ID, axis_fake_atoms[2], loc_indicator, residue_name, axis_fake_chain_ID, y_axis_ID, iCode, 0, norm, 0, occupation, B_factor, seg_ID, axis_fake_atoms[2], charge)
  z_axis <- sprintf(format, head, z_axis_ID, axis_fake_atoms[3], loc_indicator, residue_name, axis_fake_chain_ID, z_axis_ID, iCode, 0, 0, norm, occupation, B_factor, seg_ID, axis_fake_atoms[3], charge)
  
  origin <- sprintf(format, head, origin_ID, axis_fake_atoms[4], loc_indicator, residue_name, axis_fake_chain_ID, origin_ID, iCode, 0, 0, 0, occupation, B_factor, seg_ID, axis_fake_atoms[4], charge)
  
  write(x_axis, file= pdb_filepath, append= TRUE)
  write(y_axis, file= pdb_filepath, append= TRUE)
  write(z_axis, file= pdb_filepath, append= TRUE)
  write(origin, file= pdb_filepath, append= TRUE)
  
  #Connecting x, y and z to origin
  connect_x <- paste("CONECT", x_axis_ID, origin_ID)
  connect_y <- paste("CONECT", y_axis_ID, origin_ID)
  connect_z <- paste("CONECT", z_axis_ID, origin_ID)
  
  write(connect_x, file= pdb_filepath, append= TRUE)
  write(connect_y, file= pdb_filepath, append= TRUE)
  write(connect_z, file= pdb_filepath, append= TRUE)
  
  #Display CONECT elements when loading pml file
  enable_connections <- "show lines"
  write(enable_connections, file= pml_filepath, append= TRUE)
  
}
