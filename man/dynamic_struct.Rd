\name{dynamic_struct}

\Rdversion{1.1}

\alias{dynamic_struct}

\title{
  Creates the data structure for the analysis of side chain dihedral angles   
}

\description{
  Given a structure pdb file, a trajectory dcd file and frame indices, gathers information on side chain dihedral angles in a unique structure. This structure will be used in correlation/covariation methods aimed at analyzing side chain rotational motions during molecular dynamics simulations.
}

\usage{
  dynamic_struct(pdb, trj, frames)
}

\arguments{
 \item{pdb}{
   Filepath of protein description (pdb file)
 }
 \item{trj}{
   Filepath of trajectory file (dcd file)
 }
 \item{frames}{
   Stride length indicating the ratio of selected frames in dcd file
 }
}


\value{
  Returns a list of class 'structure' with the following elements containing information on the sequence, structure, trajectory and side chain dihedral angles (values and names) of the protein during the simulation:
  \item{pdb}{
  an object of class 'pdb' created by the \code{read.pdb}function from the \code{bio3d} package
  }
  \item{dcd}{
  A numeric matrix of xyz coordinates with a frame per row and a Cartesian coordinate per column. Created by the \code{read.dcd}function from the \code{bio3d} package  
  }
  \item{xyz}{
  A numeric matrix of xyz coordinates with a frame per row and a Cartesian coordinate per column. For each frame, the protein coordinates have been fitted on the pdb structure using the \code{fit.xyz} from the \code{bio3d} package  .
  }
  \item{tor}{
  A numeric matrix of side chain dihedral angles with a frame per row and a dihedral angle per column. Contains only side chain dihedral angles. Created by the \code{xyz2tor} function from the \code{bio3d} package.  
  }
  \item{tor.names}{
  a character vector with the names of all side chain chi dihedral angles. They are written as "M.chiN" where M is the residue number and N the dihedral angle chi (chi1, chi2,...). Positions for Alanine and Glycine residues which do not have side chain dihedral angles are omitted. 
  }
  \item{prot.seq}{
  a character vector with the sequence of the protein. 
  }
  \item{nb_residues}{
  a numeric value indicating the number of residues in the protein
  }
  \item{nb_frames}{
  a numeric vector indicating the frames that are used (first, last, stride)
  }



}

\author{
  Bruck TADDESE and Antoine GARNIER
}

\examples{
   #Input files
  pdb <- system.file("rotamer/toy_coordinates.pdb", package= "Bios2cor")
  trj <- system.file("rotamer/toy_dynamics.dcd", package= "Bios2cor")
  conversion_file <- system.file("rotamer/dynameomics_rotameres.csv", package= "Bios2cor")
  
  wanted_frames <- seq(from= 1, to= 40, by= 2)
  dynamic_struct <- dynamic_struct(pdb, trj, wanted_frames)
  
}

\keyword{conversion}
\keyword{angle}
\keyword{trajectory}
\keyword{protein}
 
