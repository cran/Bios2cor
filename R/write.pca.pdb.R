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
write.pca.pdb <- function(corr_pca, filepathroot=NULL, trio_comp= c(1:3)){

  if (missing(corr_pca)) {
      stop("A PCA object created by the centered_pca function is required")
  }

  pca_coord <- corr_pca$coord
  pca_1 <- trio_comp[1]
  pca_2 <- trio_comp[2]
  pca_3 <- trio_comp[3]
  
  if (is.null(filepathroot)) {
    pdb_filepath <- paste("PCA_", pca_1, "_",pca_2, "_",pca_3, ".pdb", sep ="")
    pml_filepath <- paste("PCA_", pca_1, "_",pca_2, "_",pca_3, ".pml", sep ="")
  } else {
    pdb_filepath <- paste(filepathroot, "_PCA_", pca_1, "_",pca_2, "_",pca_3, ".pdb", sep ="")
    pml_filepath <- paste(filepathroot, "_PCA_", pca_1, "_",pca_2, "_",pca_3, ".pml", sep ="")   
  }


  pca_positions <- rownames(pca_coord)
  pca_size <- length(pca_positions)


#check whether data are from a MSA or a trajectory
 numbering <- grep(".", pca_positions, fixed= TRUE)

  alpha <- c("A","B","C","D")
  colors <- c("blue", "red", "green", "orange", "white")
  background_color <- paste("bg_color", colors[5])
  write(background_color, file= pml_filepath, append= FALSE)
  for(i in 1:4){
      chain_coloration <- paste("color", colors[i], ", chain", alpha[i])
      write(chain_coloration, file= pml_filepath, append= TRUE)
  }



#Configuring elements
  as_sphere <- "as sphere"
  sphere_scale <- "set sphere_scale, 1"
  write(as_sphere, file= pml_filepath, append= TRUE)
  write(sphere_scale, file= pml_filepath, append= TRUE)




## Writing PDB file

  pdb_line <- "REMARK PDB file created for visualization of PCA analysis"
  write(pdb_line, file= pdb_filepath, append= FALSE)
 
##fields required in each line
  #fields without change
  head <- "HETATM"
  atom_name <- "O"
  res_type <- "HOH"
  occupation <- "1.00"
  B_factor <- 100.00
  atom_type <- "O"
  format <- "%6s%5s  %-4s%3s %1s%4s    %8s%8s%8s   %4s%6s     %4s"  

  #fields with changes 

  if(length(numbering) == 0){

   for(pos in 1:pca_size){
      x <- pca_coord[pos, pca_1]
      y <- pca_coord[pos, pca_2]
      z <- pca_coord[pos, pca_3]

     if (!is.na(x) & !is.na(y) & !is.na(z)) { 
      residue_number <- pca_positions[pos]
      chain_ID <- "A"
      pdb_line <- sprintf(format, head, pos, atom_name, res_type, chain_ID, residue_number, round(x,digits=3), round(y,digits=3), round(z,digits=3), occupation, B_factor, atom_type)
      write(pdb_line, file= pdb_filepath, append= TRUE)
     }
    }
   last_number <- as.numeric(residue_number)

  }else{

  resno <- sub("\\..*$", "", pca_positions)
  angle <- sub("^[0-9]*\\.","", pca_positions)

  for(pos in 1:pca_size){
      x <- pca_coord[pos, pca_1]
      y <- pca_coord[pos, pca_2]
      z <- pca_coord[pos, pca_3]

    if (!is.na(x) & !is.na(y) & !is.na(z)) { 
    residue_number <- resno[pos]
    angle_ID <- angle[pos]

    if (angle_ID == "chi1") {chain_ID <- "A"}
    if (angle_ID == "chi2") {chain_ID <- "B"}
    if (angle_ID == "chi3") {chain_ID <- "C"}
    if (angle_ID == "chi4") {chain_ID <- "D"}

    pdb_line <- sprintf(format, head, pos, atom_name, res_type, chain_ID, residue_number, round(x,digits=3), round(y,digits=3), round(z,digits=3), occupation, B_factor, atom_type)
    
    write(pdb_line, file= pdb_filepath, append= TRUE)
   }
  }
    last_number <- as.numeric(residue_number)


}
  
  #Creating elements used to create axis
  max_x <- max(abs(na.omit(pca_coord[,1])))
  norm <- round(max_x*2, digits=3)
  axis_fake_atoms <- c("O","C","N", "Pb")
  axis_fake_chain_ID <- "Z"

  x_pos <- pca_size+1
  y_pos <- pca_size+2
  z_pos <- pca_size+3
  origin_pos <- pca_size+4

  x_residue <- last_number+1
  y_residue <- last_number+2
  z_residue <- last_number+3
  origin_residue <- last_number+4
 
  x_name <- "XXX"
  y_name <- "YYY"
  z_name <- "ZZZ"
  origin_name <- "000"
  
  x_axis <- sprintf(format, head, x_residue, axis_fake_atoms[1], x_name, axis_fake_chain_ID, x_residue, norm, 0, 0, occupation, B_factor, axis_fake_atoms[1])
  y_axis <- sprintf(format, head, y_residue, axis_fake_atoms[2], y_name, axis_fake_chain_ID, y_residue, 0, norm, 0, occupation, B_factor, axis_fake_atoms[2])
  z_axis <- sprintf(format, head, z_residue, axis_fake_atoms[3], z_name, axis_fake_chain_ID, z_residue, 0, 0, norm, occupation, B_factor, axis_fake_atoms[3])
  origin <- sprintf(format, head, origin_residue, axis_fake_atoms[4], origin_name, axis_fake_chain_ID, origin_residue, 0, 0, 0, occupation, B_factor, axis_fake_atoms[4])
  
  write(x_axis, file= pdb_filepath, append= TRUE)
  write(y_axis, file= pdb_filepath, append= TRUE)
  write(z_axis, file= pdb_filepath, append= TRUE)
  write(origin, file= pdb_filepath, append= TRUE)
  
  #Connecting x, y and z to origin
  connect_x <- paste("CONECT", x_residue, origin_residue)
  connect_y <- paste("CONECT", y_residue, origin_residue)
  connect_z <- paste("CONECT", z_residue, origin_residue)
  
  write(connect_x, file= pdb_filepath, append= TRUE)
  write(connect_y, file= pdb_filepath, append= TRUE)
  write(connect_z, file= pdb_filepath, append= TRUE)
  write("END", file= pdb_filepath, append= TRUE)

  
  #Display CONECT elements when loading pml file
  enable_connections <- "show lines"
  write(enable_connections, file= pml_filepath, append= TRUE)
  hidden_spheres <- "hide spheres, chain Z" 
  write(hidden_spheres, file= pml_filepath, append= TRUE)
 
}
