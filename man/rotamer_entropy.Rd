\name{rotamer_entropy}

\Rdversion{1.1}

\alias{rotamer_entropy}

\title{
Entropy score
}

\description{
  Measures a "dynamic entropy" score of each dihedral angle based on the number of rotameric changes during the molecular dynamics trajectory.
}

\usage{
  rotamer_entropy(rotamers)
}

\arguments{
 \item{rotamers}{
  A character matrix of type 'rotamers' that is produced by the \link{angle_conversion} function. 
  The matrix indicates the rotameric state of each side chain dihedral angle for each frame of the trajectory.
 }
}

\details{
  The entropy score S is computed by summing the number of rotameric changes over all frames, normalized to the number of frames.
}

\value{
  A numeric vector containing a "dynamic entropy" score for each side chain dihedral angle during the trajectory. The score is comprised between 0 (no change in the rotameric state during the trajectory) and 1 (rotameric change for every frame of the trajectory).   
}

\author{
  Antoine GARNIER
}

\examples{
 #Reading pdb and dcd files
  pdb <- system.file("rotamer/toy_coordinates.pdb", package= "Bios2cor")
  trj <- system.file("rotamer/toy_dynamics.dcd", package= "Bios2cor")
 
  #Reading conversion file
  conversion_file <- system.file("rotamer/dynameomics_rotameres.csv", package= "Bios2cor")

  #Creating angle_conversion object
  wanted_frames <- seq(from= 1, to= 40, by= 5)
  dynamic_struct <- dynamic_struct(pdb, trj, wanted_frames)
  
  my_rotamers <- angle_conversion(dynamic_struct, conversion_file)
 
  #creating rotamer entropy object
  rotamer_entropy <- rotamer_entropy(my_rotamers)
}

\keyword{entropy}
\keyword{rotamers}
