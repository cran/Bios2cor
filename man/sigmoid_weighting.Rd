\name{sigmoid_weighting}

\Rdversion{1.1}

\alias{sigmoid_weighting}

\title{
  Sigmoid weighting for each position
}

\description{
  Given an entropy object (result of the \code{\link{entropy}} or of the \code{\link{rotamer_entropy}}function), creates a vector with a sigmoid weighting of each element based on the entropy value. The vector will be used to give different weights to each element in the \code{\link{centered_pca}} function.
}

\usage{
  sigmoid_weighting(entropy, L= 10, inf_pt= 0.1, stability= TRUE)
}

\arguments{
 \item{entropy}{
  An object created by the \code{\link{entropy}} function
 }
 \item{L}{
  A value indicating the slope on the inflexion point. A value between 5 and 10 is recommanded. Default is 10
 }
 \item{inf_pt}{
  A value between 0 and 1 corresponding to the sigmoid curve inflexion point. Default is 0.1
 }
 \item{stability}{
  Logical value indicating whether sigmoid function favors low (TRUE) or high (FALSE) entropy values. Default is TRUE
 }
}

\details{
  The object returned by the \code{\link{entropy}} function contains an entropy score for each position.
  The weighting of each position is calculated as follows :
     \deqn{ weighting[i] = \left\{ \begin{array}{ll}  \frac{1}{(1+e^(-L*( entropy[i]- inf_pt )))}, & "stability" = FALSE \\ \\  1-(\frac{1}{(1+e^(-L*( entropy[i] - inf_pt)))}), & "stability" = TRUE \end{array}  \right. }
This option is useful to limit analysis to side chain dihedral angles with slow or fast moving angles (TRUE or FALSE stability option, respectively).
}

\value{
   A numeric matrix that contains a weighting score for each position in the alignment based on sigmoid function.

}

\author{
  Antoine GARNIER
}

\examples{
  align <- import.msf(system.file("msa/toy_align.msf", package = "Bios2cor"))

  #Creating ENTROPY object
  entropy <- entropy(align)
  
  filter <- sigmoid_weighting(entropy, L= 10, inf_pt= 0.1, stability= TRUE)
}

\keyword{weighting}
\keyword{entropy}
\keyword{sigmoid}
