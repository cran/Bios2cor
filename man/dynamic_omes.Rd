\name{dynamic_omes}

\Rdversion{1.1}

\alias{dynamic_omes}

\title{
  OMES (Observed minus Expected Squared) function applied to rotamers in molecular dynamics simulations
}

\description{
  Calculates difference between observed and expected occurrences of each possible pair of rotamers (x, y) for i and j dihedral angles over all frames
}

\usage{
  dynamic_omes(
    dynamic_structure,
    rotamers,
    res_selection= 
      c("C","I","L","M","V","R","H","K","D","E","N","Q","F","Y","W","T","S","P")
  )
}

\arguments{
 \item{dynamic_structure}{
   An object of class 'structure' that is created by the \code{\link{dynamic_structure}} function
 }
 \item{rotamers}{
 A character matrix that is produced by the \link{angle2rotamer} function. The matrix indicates the rotameric state of each side chain dihedral angle for each frame of the trajectory.
 }
 \item{res_selection}{
  List of amino acids that will be taken into account in the correlation/covariation matrix. By default, all the amino acids are taken into account except Gly and Ala, with no side chain dihedral angles.  
 }
}

\details{
  The OMES score for angles [i,j] has been computed with the following formula :
  
  \deqn{{OMES(i,j)} = \frac{1}{N} \sum_{x,y}^{ }(N_{x,y}^{obs}(i,j)-N_{x,y}^{ex}(i,j))^2}
  
  with : \emph{\eqn{N_{x,y}^{ex}(i,j) = p_{x}(i)p_{y}(j)N}}

  where :
    \itemize{
      \item {\emph{\eqn{N_{x,y}^{obs}(i,j)}}}{ is number of times that each (x,y) rotamer pair is observed at angles i and j}
      \item {\emph{\eqn{N_{x,y}^{ex}(i,j)}}}{ is number of times that each (x,y) rotamer pair would be expected, based on the frequency of rotamer x and y at angles i and j}
      \item {\emph{\eqn{N}}}{ is the number of frames}
      \item {\emph{\eqn{p_{x}(i)}}}{ is the frequency of rotamer x at angle i}
      \item {\emph{\eqn{p_{y}(j)}}}{ is the frequency of rotamer y at angle j}
    }
}

\value{
  returns a list of four elements which are numeric matrices containing (1) the correlation/covariation scores for each pair of rotamers (score), (2) the Z-scores for each pair of rotamers (Zscore), (3) the correlation/covariation scores for each pair of rotamers with zero values for autocorrelation (correlation within the same side chain) (score_noauto) and (4) the Z-scores calculated without autocorrelation pairs and zero values for autocorrelation pairs (Zscore_noauto). 
}

\author{
   Antoine GARNIER, Lea BELLENGER, and Marie CHABBERT 
}

\references{
Fodor AA and Aldrich RW. Influence of conservation on calculations of amino acid covariance in multiple sequence alignments. Proteins. 2004;56(2):211-21.
}

\examples{
  #Reading pdb and dcd files
  pdb <- system.file("rotamer/toy_coordinates.pdb", package= "Bios2cor")
  trj <- system.file("rotamer/toy_dynamics.dcd", package= "Bios2cor")

  #Creating dynamic_structure object
  wanted_frames <- seq(from = 5, to = 40, by = 10)
  dynamic_structure <- dynamic_structure(pdb, trj, wanted_frames)

  #Creating rotamers object using conversion_file
  conversion_file <- system.file("rotamer/dynameomics_rotamers.csv", package= "Bios2cor")
  rotamers <- angle2rotamer(dynamic_structure, conversion_file)

  #Creating correlation object for selected residues with OMES method
  wanted_residues <- c("W")
  omes_corr <- dynamic_omes(dynamic_structure, rotamers, wanted_residues)
}


