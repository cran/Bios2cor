\name{cyto_entropy}

\Rdversion{1.1}

\alias{cyto_entropy}

\title{
  Creation of a entropy file in Cytoscape format
}

\description{
  Given an entropy object, stores entropy scores in a file formatted as a Cytoscape node attribute file
}

\usage{
  cyto_entropy(entropy, contact= NULL, filepath)
}

\arguments{
 \item{entropy}{
  An object created by an entropy function
 }
 \item{contact}{
  Result of \code{\link{corr_contact}} function that contains top positions to be taken into account
 }
 \item{filepath}{
  The full path name of the output file
 }
}

\details{
  This function purpose is storing the entropy scores in a file that will read by Cytoscape as node attributes
}

\author{
  Antoine GARNIER
}

\examples{
  align <- import.msf(system.file("msa/toy_align.msf", package = "Bios2cor"))

  #Creating ENTROPY object and writing results in cytoscape format
  entropy <- entropy(align)
  
  cyto_ent <- cyto_entropy(entropy, contact = NULL, "~/cyto_entropy.txt")
}

\keyword{cytoscape}
\keyword{entropy} 
