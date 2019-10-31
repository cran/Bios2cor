\name{dynamic_entropy}

\Rdversion{1.1}

\alias{dynamic_entropy}

\title{
Variability score
}

\description{
  Measures a "dynamic entropy" or variability score of each dihedral angle based on the number of rotameric changes during the molecular dynamics trajectory.
}

\usage{
  dynamic_entropy(rotamers)
}

\arguments{
 \item{rotamers}{
  A character matrix of type 'rotamers' that is produced by the \code{angle2rotamer} function. 
  The matrix indicates the rotameric state of each side chain dihedral angle for each frame of the trajectory.
 }
}

\details{
  The "dynamic entropy" score S is computed by summing the number of rotameric changes over all frames, normalized to the number of frames. It is not a "true entropy" score but gives usefull information on variability of the dihedral angle during the MD simulation.
}

\value{
  A numeric vector containing a "dynamic entropy" score for each side chain dihedral angle during the trajectory. The score is comprised between 0 (no change in the rotameric state during the trajectory) and 1 (rotameric change for every frame of the trajectory).   
}

\author{
 Antoine GARNIER, Lea BELLENGER and Marie CHABBERT 
}

\examples{
  #Reading pdb and dcd files
  pdb <- system.file("rotamer/toy_coordinates.pdb", package= "Bios2cor")
  trj <- system.file("rotamer/toy_dynamics.dcd", package= "Bios2cor")
 
  #Reading conversion file
  conversion_file <- system.file("rotamer/dynameomics_rotamers.csv", package= "Bios2cor")

  #Creating the dynamic_structure and rotamers objects
  wanted_frames <- seq(from = 1, t = 40, by = 5)
  dynamic_structure <- dynamic_structure(pdb, trj, wanted_frames)
  rotamers <- angle2rotamer(dynamic_structure, conversion_file)
 
  #creating the dynamic_entropy object
  dynamic_entropy <- dynamic_entropy(rotamers)
}


