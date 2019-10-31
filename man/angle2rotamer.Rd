\name{angle2rotamer}

\Rdversion{1.1}

\alias{angle2rotamer}

\title{
  Converts dihedral angles to rotamers
}

\description{
  Given an object of class 'structure' and an angle to rotamer conversion file, associates a rotamer to each dihedral angle value.
  The object of class 'structure' contains dihedral angle values for each side chain dihedral angle and each frame of the trajectory.
  The conversion file is a reference file that contains the rotamer to be associated with a dihedral angle value, depending on the residue type and the dihedral angle considered. This function will allow to compute correlation/covariation scores between rotameric states.
}

\usage{
  angle2rotamer(dynamic_structure, 
    conversion_file=system.file("rotamer/dynameomics_rotamers.csv", package= "Bios2cor"))
}

\arguments{
 \item{dynamic_structure}{
   Dihedral angle structure, result of the \code{dynamic_structure} function
 }
 \item{conversion_file}{
   The file containing the rotamer to be associated with each residue dihedral angle depending of the dihedral angle value.

   Each line contains five fields, separated by ','. The five fields represent the amino acid name ("R", "N",...), the dihedral angle name("chi1", "chi2",...), the associated rotamer ("g+", "t", "g-"), the start and stop angles (between -180 and 180). For example, for the \emph{chi1} angle of the valine residue, a torsion angle between 0 and 120 is associated to rotameric state \emph{g+}.

   Default is the "rotamer/dynameomics_rotamers.csv" conversion file provided with the Bios2cor package from the dynameomics database (\code{http://dynameomics.org}) .
 }
}

\details{
  In the torsion object and in the conversion file, dihedral angle values vary between -180 and 180.
}

\value{
  A character matrix containing the rotameric state of each side chain dihedral angle for each frame in the trajectory, depending on the dihedral angle value and on a conversion file.
}

\author{
  Antoine GARNIER and Lea BELLENGER
}

\references{
 Van der Kamp MW, Schaeffer RD, Jonsson AL, Scouras AD, Simms AM, Toofanny RD, Benson NC, Anderson PC, Merkley ED, Rysavy S, Bromley D, Beck DAC, and Daggett V. Dynameomics: A comprehensive database of protein dynamics. Structure, 18: 423-435, 2010. 
}


\examples{
  #Reading pdb and dcd files
  pdb <- system.file("rotamer/toy_coordinates.pdb", package= "Bios2cor")
  trj <- system.file("rotamer/toy_dynamics.dcd", package= "Bios2cor")
 
  #Reading conversion file 
  #conversion_file <- system.file("rotamer/dynameomics_rotamers.csv", package= "Bios2cor")

  #Creating dynamic_structure object
  wanted_frames <- seq(from = 1, to = 40, by = 2)
  dynamic_structure <- dynamic_structure(pdb, trj, wanted_frames)
  
  #Creating rotamers object with default conversion_file
  rotamers <- angle2rotamer(dynamic_structure)
  
  #Creating rotamers object with conversion_file   
  #rotamers <- angle2rotamer(dynamic_structure, conversion_file)
 
}


