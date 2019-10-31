\name{dynamic_structure}

\Rdversion{1.1}

\alias{dynamic_structure}

\title{
  Creates the data structure for the analysis of side chain dihedral angles   
}

\description{
  Given a structure pdb file, a trajectory dcd file and frame indices, gathers information on side chain dihedral angles in a unique structure. This structure will be used in correlation/covariation methods aimed at analyzing side chain rotational motions during molecular dynamics simulations.
}

\usage{
  dynamic_structure(pdb, trj, frames=NULL)
}

\arguments{
 \item{pdb}{
   Filepath of the  pdb file
 }
 \item{trj}{
   Filepath of trajectory file (dcd file)
 }
 \item{frames}{
   Indicates the selected frames for the analysis, created with the \code{seq} function (Usage: frames <-seq(from ,to , by= ). Default is NULL (all the frames of the trajectory are taken into account).  
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
  A numeric matrix of xyz coordinates with a frame per row and a Cartesian coordinate per column. For each frame, the protein coordinates have been fitted on the pdb structure using the \code{fit.xyz} from the \code{bio3d} package  
  }
  \item{tor}{
  A numeric matrix of side chain dihedral angles with a frame per row and a dihedral angle per column. Contains only side chain dihedral angles. Created by the \code{xyz2tor} function from the \code{bio3d} package  
  }
  \item{tor.names}{
  a character vector with the names of all side chain chi dihedral angles. They are written as "M.chiN" where M is the residue number and N the dihedral angle chi (chi1, chi2,...). Alanine and Glycine residues which do not have side chain dihedral angles are omitted. 
  }
  \item{tor.resno}{
  a character vector with the residue number M of all side chain chi dihedral angles. 
  }
  \item{tor.angle}{
  a character vector with the dihedral angle (chiN) of all side chain chi dihedral angles. 
  }
  \item{nb_torsions}{
  a numeric value indicating the total number of dihedral angles 
  }
  \item{prot.seq}{
  a character vector with the sequence of the protein 
  }
  \item{nb_residues}{
  a numeric value indicating the number of residues in the protein
  }
  \item{nb_frames}{
  a numeric value indicating the total number of selected frames 
  }
  \item{frames}{
  a numeric vector indicating the selected frames 
  }



}

\author{
  Bruck TADDESE, Antoine GARNIER and Marie CHABBERT
}

\examples{
  #Reading pdb and dcd files
  pdb <- system.file("rotamer/toy_coordinates.pdb", package= "Bios2cor")
  trj <- system.file("rotamer/toy_dynamics.dcd", package= "Bios2cor")
  
  #Creating dynamic_structure object for selected frames  
  wanted_frames <- seq(from = 1, to = 40, by = 2)
  dynamic_structure <- dynamic_structure(pdb, trj, wanted_frames)

  #Creating dynamic_structure object for all the frames  
  #dynamic_structure <- dynamic_structure(pdb, trj)
}

 
