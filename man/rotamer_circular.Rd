\name{rotamer_circular}

\Rdversion{1.1}

\alias{rotamer_circular}

\title{
  Circular correlation
}

\description{
  Calculates circular correlation scores between side chain dihedral angles during a molecular dynamics trajectory. 
}

\usage{
  rotamer_circular(
    dynamic_struct,
    res_selection=
      c("C","I","L","M","V","R","H","K","D","E","N","Q","F","Y","W","T","S","P"),
    z_score= TRUE,
    auto_pairing= FALSE
  )
}

\arguments{
 \item{dynamic_struct}{
   Object of class 'structure' that is created by the \bold{\emph{dynamic_struct}} function
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
  This function uses the \emph{cor.circular} function from the \emph{circular} package based on a circular version of the Pearson coeefficient. 
}

\value{
  A list of two elements : a matrix containing a covariation score for each pair of rotamer and its normalized version
}

\author{
   Bruck TADESSE and Antoine GARNIER
}

\references{
  Circular Statistics, from ``Topics in circular Statistics'' (2001) S. Rao Jammalamadaka and A. SenGupta, World Scientific.
}


\examples{
  #Reading pdb file
  pdb <- system.file("rotamer/toy_coordinates.pdb", package= "Bios2cor")
  trj <- system.file("rotamer/toy_dynamics.dcd", package= "Bios2cor")

  #conversion_file <- system.file("rotamer/dynameomics_rotameres.csv", package= "Bios2cor")

  wanted_residues <- c("H","N","Q","F","Y","W")
  
  wanted_frames <- seq(from= 1, to= 40, by= 2)
  dynamic_struct <- dynamic_struct(pdb, trj, wanted_frames)
    
  dihed_corr <- rotamer_circular(dynamic_struct, wanted_residues)
}

\keyword{circular}
\keyword{correlation}
\keyword{protein}
