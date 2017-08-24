\name{create_entropyfile}

\Rdversion{1.1}

\alias{create_entropyfile}

\title{
  Stores entropy values
}

\description{
  Given an entropy object, stores each element (position or dihedral angle) and its entropy value in a csv file.
}

\usage{
  create_entropyfile(entropy, filepath)
}

\arguments{
 \item{entropy}{
  An object created by the \bold{\emph{entropy}} or the \code{\link{rotamer_entropy}}function which required an alignment file name
 }
 \item{filepath}{
  The full path name of the output file 
 }
}

\details{
  The \bold{\emph{entropy}} function calculates entropy score for each position of an alignment file. The \code{\link{rotamer_entropy}}function calculate a "dynamic entropy" score for each side chain dihedral angle of a protein during a molecular simulations. 
}


\value{
  A csv file containing entropy scores
}

\author{
  Antoine GARNIER
}

\examples{
   align <- import.msf(system.file("msa/toy_align.msf", package = "Bios2cor"))

  #Creating ENTROPY object (as an example)
  entropy <- entropy(align)
  
  create_entropyfile(entropy, "results_entropy.txt")
}

\keyword{entropy}
 
