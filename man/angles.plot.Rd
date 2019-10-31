\name{angles.plot}

\Rdversion{1.1}

\alias{angles.plot}

\title{
  plots the time evolution of the dihedral angles in top pairs
}

\description{
  Given an object of class 'structure' and names of dihedral angles, creates plots of the dihedral angles as a function of frames in a pdf file. 
}

\usage{
  angles.plot(dynamic_structure, angles, filepathroot=NULL)
}

\arguments{
 \item{dynamic_structure}{
   Dihedral angle structure, result of the \code{dynamic_structure} function
 }
 \item{angles}{
   A vector containing the names of the dihedral angles to be visualized.Default is NULL (all the torsional angles of the dynamic_structure object are taken into account). 
 }
 \item{filepathroot}{
   Root of the full path name for the output file. Default is NULL. In this case, the output file is named "ANGLES.pdf". When filepathroot is not NULL, a "_ANGLES.pdf" extention is added to the filepathroot name.
 }
}

\details{
  The object of class 'structure' contains the side chain dihedral angles (between -180 and 180) for each residue in the protein, for each frame of the molecular simulations. This function allows visualisation of the evolution of selected angles. 
}

\value{
  returns a pdf file containing the plots of the frame dependance of each element included in argument \code{angles}.
}

\author{
  Antoine GARNIER and Marie CHABBERT
}

\examples{  
  #Indicating file path for output files
  wd <- tempdir()
  #wd <- getwd()
  file <- file.path(wd,"test_dyn1") 

  #Reading pdb and dcd files
  pdb <- system.file("rotamer/toy_coordinates.pdb", package= "Bios2cor")
  trj <- system.file("rotamer/toy_dynamics.dcd", package= "Bios2cor")
  
  #Creating dynamic_structure object for wanted frames
  wanted_frames <- seq(from= 1, to= 40, by= 2)
  dynamic_structure <- dynamic_structure(pdb, trj, wanted_frames)

  #Calculating circular correlation between dihedral angles of selected residues 
  wanted_residues <- c("H","N","Q","F","Y","W")
  dihed_corr <- dynamic_circular(dynamic_structure, wanted_residues)

  #Selecting a correlation matrix
  dihed_corr <- dihed_corr$Zscore
  
  #Selecting angles of interest (here from the "top_pairs_analysis" function)
  top_angles <- top_pairs_analysis(dihed_corr, top= 25, file)
  my_angles <- unlist(top_angles$positions)
  
  #Creating plots of the time evolution of the dihedral angles 
  evol_angles <- angles.plot(dynamic_structure, my_angles, file)

}

 
