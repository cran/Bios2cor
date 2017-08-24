\name{rotamer_omes}

\Rdversion{1.1}

\alias{rotamer_omes}

\title{
  OMES(Observed minus Expected Squared) function applied to rotamers in molecular dynamics simulations
}

\description{
  Calculates difference between observed and expected occurrences of each possible pair of rotamers (x, y) for i and j dihedral angles over all frames
}

\usage{
  rotamer_omes(
    dynamic_struct,
    rotamers,
    res_selection= 
      c("C","I","L","M","V","R","H","K","D","E","N","Q","F","Y","W","T","S","P"),
    z_score= TRUE,
    auto_pairing= FALSE
  )
}

\arguments{
 \item{dynamic_struct}{
   An object of class 'structure' that is created by the \bold{\emph{dynamic_struct}} function
 }
 \item{rotamers}{
 A character matrix that is produced by the \link{angle_conversion} function. The matrix indicates the rotameric state of each side chain dihedral angle for each frame of the trajectory.
 }
 \item{res_selection}{
  Selection of amino acid types that will be taken into account in the covariation matrix. Allows to limit the analysis to a limited selection of amino acid types. 
 }
 \item{z_score}{
    A logical value for Z-score normalisation of the covariation matrix. Default is TRUE.
 }
 \item{auto_pairing}{
   A logical value that maintains (TRUE) or removes (FALSE) covariation scores between dihedral angles within a same residue in the covariation matrix. DEFAULT is FALSE.
 }
}

\details{
  The OMES score at angles [i,j] has been computed with the following formula :
  
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
  returns a list of one or two elements which are a numeric matrix containing the correlation scores for each pair of rotamers and, optionaly, a  numeric matrix containing the Z-scores 
}

\author{
   Antoine GARNIER, Madeline DENIAUD and Lea Bellenger
}

\references{
Fodor AA and Aldrich RW. Influence of conservation on calculations of amino acid covariance in multiple sequence alignments. Proteins. 2004;56(2):211-21.
}

\examples{
  #Calculating rotamers
  pdb <- system.file("rotamer/toy_coordinates.pdb", package= "Bios2cor")
  trj <- system.file("rotamer/toy_dynamics.dcd", package= "Bios2cor")

  conversion_file <- system.file("rotamer/dynameomics_rotameres.csv", package= "Bios2cor")

  wanted_residues <- c("W")
  
  nb_frames_wanted <- 40
  wanted_frames <- seq(from= 5, to= nb_frames_wanted, by= 10)
  dynamic_struct <- dynamic_struct(pdb, trj, wanted_frames)
  
  my_angles <- angle_conversion(dynamic_struct, conversion_file)
  
  #Creating OMES object
  omes_corr <- rotamer_omes(dynamic_struct, my_angles, wanted_residues, z_score= FALSE)
}

\keyword{omes}
\keyword{rotamer}
