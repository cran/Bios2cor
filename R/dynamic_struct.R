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

dynamic_struct <- function(pdb, trj, frames){
  ##Calculating rotamer informations
  #Reading pdb file
  pdb <- read.pdb(pdb)
  trj <- read.dcd(trj)

  #Selecting only "CA" atoms
  ca.inds <- atom.select(pdb, elety = "CA")

  #Getting xyz coordinates using fit.xyz form bio3d package
  xyz <- fit.xyz(fixed = pdb$xyz, mobile = trj, fixed.inds = ca.inds$xyz,  mobile.inds = ca.inds$xyz)
  nb_frames <- length(xyz[,1])
  
  #Creating torsion object using xyz2torsion function from bio3d package
  if(is.null(frames)){
    tor <- xyz2torsion(pdb, xyz, tbl = "sidechain", ncore= 1)
  } else {
    tor <- xyz2torsion(pdb, xyz[frames,], tbl = "sidechain", ncore= 1)
  }
  nb_residus <- length(tor[1,])
  nb_frames <- length(tor[,1])
  
  #Associating special letter for each residue, using pdbseq function
  prot.seq <- pdbseq(pdb)
  
  tor.names <- colnames(tor)
  tor.resno <- sub("\\..*$", "", tor.names)
  tor.angle <- sub("^[0-9]*\\.","", tor.names)
  tor.seq <- prot.seq[tor.resno]

  ##Filling dihedral angle structure
  res <- list()
  res$pdb <- pdb
  res$trj <- trj
  
  res$ca.inds <- ca.inds
  
  res$xyz <- xyz
  res$nb_frames <- nb_frames
  
  res$tor <- tor
  res$nb_residus <- nb_residus
  res$nb_frames <- nb_frames
  
  res$prot.seq <- prot.seq
  
  res$tor.names <- tor.names
  res$tor.resno <- tor.resno
  res$tor.angle <- tor.angle
  res$tor.seq <- tor.seq

  class (res) <- c("structure")  

  return (res)
}
