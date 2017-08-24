\name{rotamer_mip}

\Rdversion{1.1}

\alias{rotamer_mip}

\title{
  MIP(Mutual Information Product) function applied to rotamers in molecular dynamics simulations
}

\description{
  Calculates a mutual information score (MI) based on the probability of joint occurrence of events and corrects it with the average product  which is subtracted from the MI score. 
}

\usage{
  rotamer_mip(
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
  A character matrix of type 'rotamers' that is produced by the \link{angle_conversion} function. 
  The matrix indicates the rotameric state of each side chain dihedral angle for each frame of the trajectory.
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
  The MIp score at position [i,j] has been computed with the following formula :
  
  \deqn{{MIp(i,j)} = MI(i,j) - \frac{MI(i,\bar{j})MI(\bar{i},j)}{<MI>}}
  
  with :
    \itemize{
      \item {\emph{\eqn{{MI(i,j) = \sum_{x,y}^{ } p_{x,y}(i,j) ln\frac{p_{x,y}(i,j)}{p_{x}(i)p_{y}(j)}}}}} {}
      \item {\emph{\eqn{MI(i,\bar{j}) = \frac{1}{n-1} \sum_{j \neq i}^{ } MI(i,j)}}} {}
      \item {\emph{\eqn{MI(\bar{i},j) = \frac{1}{n-1} \sum_{i \neq j}^{ } MI(i,j)}}} {}
      \item {\emph{\eqn{<MI> = \frac{2}{n(n-1)} \sum_{i,j}^{ }MI(i,j)}}} {}
    }
    
  and where \emph{\eqn{p_{x,y}(i,j)}} is the frequency of the rotamer pair (x,y) at dihedral angles i and j.
      
  N.B. this formula has been widely applied in the field of sequence covariation but favors pairs with high entropy.
}

\value{
  A list of two elements : a matrix containing a correlation score for each pair of rotamer and its normalized version
}

\author{
   Antoine GARNIER and Madeline DENIAUD
}

\references{
 Dunn SD, Wahl LM, Gloor GB. Mutual information without the influence of phylogeny or entropy dramatically improves residue contact prediction. Bioinfor;atics 2008;24:333-340.
 Martin LC, Gloor GB, Dunn SD, Wahl LM. Using infor;ation theory to search for co-evolving residues in proteins. Bioinformatics 2005;21:4116-4124.
}

\examples{
  #Calculating rotamers
  pdb <- system.file("rotamer/tiny_toy_coordinates.pdb", package= "Bios2cor")
  trj <- system.file("rotamer/tiny_toy_dynamics.dcd", package= "Bios2cor")

  conversion_file <- system.file("rotamer/dynameomics_rotameres.csv", package= "Bios2cor")

  wanted_residues <- c("H","N")
  
  wanted_frames <- seq(from= 5, to= 40, by= 15)
  dynamic_struct <- dynamic_struct(pdb, trj, wanted_frames)
  
  my_angles <- angle_conversion(dynamic_struct, conversion_file)
  
  #Creating MIP object
  mip_corr <- rotamer_mip(dynamic_struct, my_angles, wanted_residues)
}

\keyword{mip}
\keyword{rotamer}
