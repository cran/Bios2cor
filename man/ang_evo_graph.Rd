\name{ang_evo_graph}

\Rdversion{1.1}

\alias{ang_evo_graph}

\title{
  Angle evolution graphs
}

\description{
  Given an object of class 'structure' and names of dihedral angles, creates graphs to monitor the evolution of the dihedral angles along the frames and organizes them in a pdf file.
}

\usage{
  ang_evo_graph(dynamic_struct, top_positions, filepath)
}

\arguments{
 \item{dynamic_struct}{
   Dihedral angle structure, result of the \bold{\emph{dynamic_struct}} function
 }
 \item{top_positions}{
   A vector containing the names of the dihedral angles in the top X pairs; Result of \emph{corr_contact} function
 }
 \item{filepath}{
   Filepath for the output file
 }
}

\details{
  The object of class 'structure' contains the side chain dihedral angles (between -180 and 180) for each residue in the protein, for each frame of the molecular simulations. This function allows visualisation of the evolution of selected angles 
}

\value{
  returns a pdf file containing the graphs of the frame dependance of each element included in argument \code{top_positions} 
}

\author{
  Antoine GARNIER
}

\examples{  
  #Reading pdb file
  pdb <- system.file("rotamer/toy_coordinates.pdb", package= "Bios2cor")
  trj <- system.file("rotamer/toy_dynamics.dcd", package= "Bios2cor")
  
   
  #output file
  filepath <- "angles.pdf"
  
  #Calculating dynamic data (dihedral angles)
  wanted_frames <- seq(from= 1, to= 40, by= 2)
  dynamic_struct <- dynamic_struct(pdb, trj, wanted_frames)

  wanted_residues <- c("H","N","Q","F","Y","W")
  
    
  dihed_corr <- rotamer_circular(dynamic_struct, wanted_residues)
  dihed_corr <- dihed_corr$normalized
  
  #Select positions of interest (calculated from the "corr_contact" function or list)
  top_angles <- corr_contact(dihed_corr, top= 25, "contact_file.txt","top_scores.txt")
  my_angles <- unlist(top_angles$positions)
  
  #Creating ang_evol graph
  evol_angles <- ang_evo_graph(dynamic_struct, my_angles, filepath)
}

\keyword{graph}
\keyword{top}
\keyword{angle}
\keyword{trajectory}
 
