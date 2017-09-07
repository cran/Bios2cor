\name{cyto_zscore}

\Rdversion{1.1}

\alias{cyto_zscore}

\title{
  Creation of a Z-score file in Cytoscape format
}

\description{
  Given a correlation object, stores Z-scores in a file formatted as a Cytoscape edge attribute file
}

\usage{
  cyto_zscore(corr_struct, contact= NULL, filepath)
}

\arguments{
 \item{corr_struct}{
  An object created by a correlation function (\code{\link{omes}}, \code{\link{mip}}, \code{\link{elsc}}, \code{\link{mcbasc}}, \code{\link{rotamer_circular}}, \code{\link{rotamer_omes}}, \code{\link{rotamer_mip}})
 }
 \item{contact}{
  Result of \code{\link{corr_contact}} function that contains top positions to be taken into account
 }
 \item{filepath}{
  The path where the file will be created
 }
}

\details{
  This function purpose is storing Z-scores in a file that will read by Cytoscape as node attributes
}

\author{
  Antoine GARNIER
}

\examples{
   align <- import.msf(system.file("msa/toy_align.msf", package = "Bios2cor"))

  #Creating OMES object and writing Z-scores in cytoscape format
  omes <- omes(align)
  
  results <- cyto_zscore(omes, contact = NULL, "cyto_Zscore.txt")
}

\keyword{cytoscape}
\keyword{correlation}
