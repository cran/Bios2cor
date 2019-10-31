\name{xyz2torsion}
\alias{xyz2torsion}
\title{ Convert Cartesian Coordinates to Torsion Angles }
\description{
  Convert cartesian coordinate matrix to torsion angles with 
  function \code{\link{torsion.pdb}}.
}
\usage{
xyz2torsion(pdb, xyz, tbl = c("basic", "mainchain", 
     "sidechain", "all", "phi", "psi", paste("chi", 1:5, sep="")), ncore = NULL)
}
\arguments{
  \item{pdb}{ A PDB structure object as obtained from \code{\link{read.pdb}}. }
  \item{xyz}{ Cartesian coordinates as a Mx(3N) matrix. }
  \item{tbl}{ Names of torsion angles to calculate. }
  \item{ncore}{ Number of CPU cores used to do the calculation. By default (NULL), use
    all detected CPU cores. }
}
\details{
  Available values for the argument \sQuote{tbl} include:
  \itemize{
     \item Basic: "phi", "psi", "chi1".
     \item Mainchain: "phi", "psi".
     \item Sidechain: "chi1", "chi2", "chi3", "chi4", "chi5".
     \item All: all of the above.
     \item Each individual angle.
  }
}
\value{
  Returns a MxP matrix, where M is the number of frames and P the number of valid torsion angles.
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.
}
\author{ Xin-Qiu Yao }
\seealso{ \code{\link{torsion.xyz}}, \code{\link{torsion.pdb}}
 }
\note{
 New function from the bio3d package, available at <https://github.com/Grantlab/bio3d/blob/master/new_funs/xyz2torsion.R>
}
\examples{
  pdb <- system.file("rotamer/toy_coordinates.pdb", package= "Bios2cor")
  trj <- system.file("rotamer/toy_dynamics.dcd", package= "Bios2cor")
  
  pdb <- read.pdb(pdb)
  trj <- read.dcd(trj)
 
  #Selecting only "CA" atoms
  ca.inds <- atom.select(pdb, elety = "CA")

  #Getting xyz coordinates using fit.xyz form bio3d package
  xyz <- fit.xyz(fixed = pdb$xyz, mobile = trj,fixed.inds=ca.inds$xyz,mobile.inds=ca.inds$xyz)
 
  frames <- seq(from= 1, to= 40, by= 2)
  
  #Creating torsion object for side chains using xyz2torsion function from bio3d package
  tor <- xyz2torsion(pdb, xyz[frames,], tbl = "sidechain", ncore= 1)
}


