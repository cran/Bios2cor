\name{write_pdb}

\Rdversion{1.1}

\alias{write_pdb}

\title{
  PDB and PML file creation for 3D representation of PCA analysis
}

\description{
  Given a PCA structure, creates .pdb and .pml files for 3D visualization with Pymol
}

\usage{
  write_pdb(pca, trio_comp= c(1:3), pdb_filepath, pml_filepath)
}

\arguments{
 \item{pca}{
   An object created by the \code{\link{centered_pca}} function
 }
 \item{trio_comp}{
  A numeric vector of length 3 indicating the principal components to be displayed. Default is c(1, 2, 3). 
 }
 \item{pdb_filepath}{
   A string of characters indicating the full path name of the output PDB file. 
 }
 \item{pml_filepath}{
    a string of characters indicating the full path name of the output PML file for visualization with Pymol. 
 }
}

\details{
  This function creates PDB and PML files to visualize the positions of the elements (sequence positions or dihedral angles) in a 3D space corresponding to three selected components of the PCA space. The PDB file can be viewed in any molecular graphics softaware. The PML file allows a nice representation with Pymol (axis, background color, size of points and for GPCRs, color code for helices).
}

\value{
  Returns two files:  a PDB file which contains three PCA coordinates for each element in PDB format and a PML file for nice visualization with Pymol. 
}

\author{
  Antoine GARNIER
}

\examples{
  align <- import.msf(system.file("msa/toy_align.msf", package = "Bios2cor"))

  #Creating OMES object
  omes <- omes(align)
  omes <- omes$normalized
  
  # Creating PCA structures for OMES method 
  pca <- centered_pca(omes, m=NULL, pc= NULL, dec_val= 5,eigenvalues_csv= NULL)

  #Creating PDB and PML files (open PDB file with Pymol then "File > Run" PML file)
  indices <- c(1,2,3)
  write_pdb(pca, indices, "pdb_file.pdb", "pml_file.pml")
}

\keyword{pdb}
\keyword{pml}
\keyword{Pymol} 
