\name{write.pca.pdb}

\Rdversion{1.1}

\alias{write.pca.pdb}

\title{
  PDB and PML file creation for 3D representation of PCA analysis
}

\description{
  Given a PCA structure, creates .pdb and .pml files for 3D visualization with Pymol
}

\usage{
  write.pca.pdb(corr_pca, filepathroot=NULL, trio_comp= c(1:3))
}

\arguments{
 \item{corr_pca}{
   An object created by the \code{\link{centered_pca}} function from a correlation/covariation matrix
 }
 \item{filepathroot}{
   The root for the full path name of the output PDB and PML files. Default is NULL (Two PCA_l_m_n.pdb and PCA_l_m_n.pml files are created where l, m, n are the 3 selected components). If not null, '_PCA_l_m_n.pdb' and '_PCA_l_m_n.pml' extensions are added to the root name.   
 }
 \item{trio_comp}{
  A numeric vector of length 3 indicating the principal components to be displayed. Default is c(1, 2, 3). 
 }
}

\details{
  This function creates PDB and PML files to visualize the positions of the elements (sequence positions or dihedral angles) in a 3D space corresponding to three selected components of the PCA space. The PDB file can be viewed in any molecular graphics softaware. The PML file allows a nice representation with Pymol (axis, background color, size of points and for GPCRs, color code for helices).
}

\value{
  Returns two files:  a PDB file which contains three PCA coordinates for each element in PDB format and a PML file for nice visualization with Pymol. 
}

\author{
  Antoine GARNIER and Marie CHABBERT
}

\examples{

# Example for MSA
  #File path for output files
  wd <- tempdir()
  #wd <-getwd() 
  file <- file.path(wd,"test_seq5") 

  #Importing MSA file
  align <- import.msf(system.file("msa/toy_align.msf", package = "Bios2cor"))

  #Creating OMES correlation object and selecting correlation matrix
  omes <- omes(align, gap_ratio = 0.2)
  cor_omes <- omes$Zscore
 
  #Creating PCA object for selected  matrix 
  pca <- centered_pca(cor_omes, filepathroot = file, filter = NULL, pc = NULL, dec_val = 5)

  #Creating PDB and PML files (open PDB file with Pymol then "File > Run" PML file)
  indices <- c(1,2,3)
  write.pca.pdb(pca, file, indices)


 
### Example for MD
  #File path for output files
  wd <- tempdir()
  #wd <-getwd() 
  file <- file.path(wd,"test_dyn5") 

  #Redaing pdb and dcd files
  pdb <- system.file("rotamer/toy_coordinates.pdb", package= "Bios2cor")
  trj <- system.file("rotamer/toy_dynamics.dcd", package= "Bios2cor")

  #Creating dynamic_structure object for selected frames
  nb_frames_wanted <- 40
  wanted_frames <- seq(from = 5, to= nb_frames_wanted, by = 10)
  dynamic_structure <- dynamic_structure(pdb, trj, wanted_frames)

  #Creating rotamers object
  conversion_file <- system.file("rotamer/dynameomics_rotamers.csv", package= "Bios2cor")
  rotamers <- angle2rotamer(dynamic_structure, conversion_file)
  
  #Creating OMES correlation object and selecting correlation matrix
  wanted_residues <- c("W","Y","F","N")
  omes <- dynamic_omes(dynamic_structure, rotamers, wanted_residues)
  cor_omes <- omes$Zscore_noauto

  #Creating PCA object for selected matrix 
  pca <- centered_pca(cor_omes, file, filter = NULL, pc = NULL, dec_val = 5)

  #Creating PDB and PML files (open PDB file with Pymol then "File > Run" PML file)
  indices <- c(1,2,3)
  write.pca.pdb(pca, file, indices)

}


