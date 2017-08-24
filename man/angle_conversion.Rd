\name{angle_conversion}

\Rdversion{1.1}

\alias{angle_conversion}

\title{
  Conversion of dihedral angles to rotamers
}

\description{
  Given an object of class 'structure' and an angle conversion file, associates a rotamer to each dihedral angle value.
  The object of class 'structure' contains dihedral angle values for each side chain dihedral angle and each frame of the trajectory.
  The conversion file is a reference file that contains the rotamer to be associated to a dihedral angle value, depending on the residue type and the dihedral angle considered.
  This function will allow to easily compare rotameric changes that occur during the simulations. 
}

\usage{
  angle_conversion(dynamic_struct, conversion_file)
}

\arguments{
 \item{dynamic_struct}{
   Rotamer structure, result of the \bold{\emph{dynamic_struct}} function
 }
 \item{conversion_file}{
   The file containing a value (rotamer) to be associated to each residue dihedral angle depending of the dihedral angle value.
   For example, for the \emph{chi1} angle of the Valine residue, a torsion angle between 0 and 120 is associated to rotameric state \emph{g+}.
   Each line contains five fields, separated by ','. The five fields represent the residue name ("R", "N",...), the dihedral angle name("chi1", "chi2",...), the associated rotamer ("g+", "t", "g-"), the start and stop angles (between -180 and 180).
 }
}

\details{
  In the torsion object and in the conversion file, dihedral angle values are between -180 and 180.
}

\value{
  A character matrix containing the rotameric state of each side chain dihedral angle for each frame in the trajectory, depending on the dihedral angle value.
}

\author{
  Antoine GARNIER and Lea BELLENGER
}

\examples{
  #Reading pdb and dcd files
  pdb <- system.file("rotamer/toy_coordinates.pdb", package= "Bios2cor")
  trj <- system.file("rotamer/toy_dynamics.dcd", package= "Bios2cor")
 
  #Reading conversion file
  conversion_file <- system.file("rotamer/dynameomics_rotameres.csv", package= "Bios2cor")

  #Creating angle_conversion object
  wanted_frames <- seq(from= 1, to= 40, by= 2)
  dynamic_struct <- dynamic_struct(pdb, trj, wanted_frames)
  
  my_rotamers <- angle_conversion(dynamic_struct, conversion_file)
  
}

\keyword{conversion}
\keyword{angle}
\keyword{trajectory}
\keyword{protein}
