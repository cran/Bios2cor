\name{dynamic_circular}

\Rdversion{1.1}

\alias{dynamic_circular}

\title{
  Circular correlation
}

\description{
  Calculates circular correlation/covariation scores between side chain dihedral angles during a molecular dynamics trajectory. 
}

\usage{
  dynamic_circular(
    dynamic_structure,
    res_selection=
      c("C","I","L","M","V","R","H","K","D","E","N","Q","F","Y","W","T","S","P")
  )
}

\arguments{
 \item{dynamic_structure}{
   Object of class 'structure' that is created by the \code{\link{dynamic_structure}} function
 }
 \item{res_selection}{
   List of amino acids that will be taken into account in the correlation/covariation matrix. By default, all the amino acids are taken into account except Gly and Ala, with no side chain dihedral angles.  
 }
}

\details{
  This function uses the \code{cor.circular} function from the \code{circular} package based on a circular version of the Pearson coeefficient. 
}

\value{
  returns a list of four elements which are numeric matrices containing (1) the correlation/covariation scores for each pair of rotamers (score), (2) the Z-scores for each pair of rotamers (Zscore), (3) the correlation/covariation scores for each pair of rotamers with zero values for autocorrelation (correlation within the same side chain) (score_noauto) and (4) the Z-scores calculated without autocorrelation pairs and zero values for autocorrelation pairs (Zscore_noauto). 
}

\author{
   Bruck TADESSE, Antoine GARNIER, and Marie CHABBERT 
}

\references{
  Circular Statistics, from ``Topics in circular Statistics'' (2001) S. Rao Jammalamadaka and A. SenGupta, World Scientific.
}


\examples{
  #Reading pdb and dcd files
  pdb <- system.file("rotamer/toy_coordinates.pdb", package= "Bios2cor")
  trj <- system.file("rotamer/toy_dynamics.dcd", package= "Bios2cor")

  #Creating dynamic_structure object for selected frames
  wanted_frames <- seq(from = 1, to = 40, by = 2)
  dynamic_structure <- dynamic_structure(pdb, trj, wanted_frames)

  #Computing circular correlation between dihedral angles of selected residues
  res_selection <- c("H","N","Q","F","Y","W")    
  dihed_corr <- dynamic_circular(dynamic_structure, res_selection)

  #Computing correlation between all dihedral angles 
  #dihed_corr <- dynamic_circular(dynamic_structure)

}


