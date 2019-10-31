\name{dynamic_mi}

\Rdversion{1.1}

\alias{dynamic_mi}

\title{
  Mutual Information (MI) function applied to rotamers in molecular dynamics simulations
}

\description{
  Calculates a mutual information score (MI) based on the probability of joint occurrence of events.
}

\usage{
  dynamic_mi(
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
  A character matrix of type 'rotamers' that is produced by the \code{\link{angle2rotamer}} function. 
  The matrix indicates the rotameric state of each side chain dihedral angle for each frame of the trajectory.
 }
 \item{res_selection}{
   List of amino acids that will be taken into account in the correlation/covariation matrix. By default, all the amino acids are taken into account except Gly and Ala, with no side chain dihedral angles.  
 }
}

\details{
  The MI score at position [i,j] has been computed with the following formula :
  
    \itemize{
      \item {\emph{\eqn{{MI(i,j) = \sum_{x,y}^{ } p_{x,y}(i,j) ln\frac{p_{x,y}(i,j)}{p_{x}(i)p_{y}(j)}}}}} {}
     }
    
 where \emph{\eqn{p_{x,y}(i,j)}} is the frequency of the rotamer pair (x,y) at dihedral angles i and j.
      
}

\value{
 returns a list of four elements which are numeric matrices containing (1) the correlation/covariation scores for each pair of rotamers (score), (2) the Z-scores for each pair of rotamers (Zscore), (3) the correlation/covariation scores for each pair of rotamers with zero values for autocorrelation (correlation within the same side chain) (score_noauto) and (4) the Z-scores calculated without autocorrelation pairs and zero values for autocorrelation pairs (Zscore_noauto). 
}

\author{
   Antoine GARNIER and Marie CHABBERT
}

\references{
 Dunn SD, Wahl LM, Gloor GB. Mutual information without the influence of phylogeny or entropy dramatically improves residue contact prediction. Bioinfor;atics 2008;24:333-340.
 Martin LC, Gloor GB, Dunn SD, Wahl LM. Using infor;ation theory to search for co-evolving residues in proteins. Bioinformatics 2005;21:4116-4124.
}

\examples{
  #Reading pdb and dcd files
  pdb <- system.file("rotamer/tiny_toy_coordinates.pdb", package= "Bios2cor")
  trj <- system.file("rotamer/tiny_toy_dynamics.dcd", package= "Bios2cor")

  #Creating dynamic_structure object
  wanted_frames <- seq(from = 5, to = 40, by = 15)
  dynamic_structure <- dynamic_structure(pdb, trj, wanted_frames)

  #Creating rotamers object using conversion_file
  conversion_file <- system.file("rotamer/dynameomics_rotamers.csv", package= "Bios2cor")
  rotamers <- angle2rotamer(dynamic_structure, conversion_file)

  #Creating correlation object for selected residues using MI method
  wanted_residues <- c("H","N")
  mi_corr <- dynamic_mi(dynamic_structure, rotamers, wanted_residues)
}


