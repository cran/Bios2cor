\name{delta_filter}

\Rdversion{1.1}

\alias{delta_filter}

\title{
  Creation of an entropy/dynamic score delta filter for each element
}

\description{
  Given an entropy object (result of the \code{\link{entropy}} or of the \code{\link{dynamic_entropy}} function), creates a vector with a delta filter of each element based on the entropy value. The vector will be used to limit the analysis to the elements in the given entropy range in the \code{\link{centered_pca}} and \code{\link{top_pairs_analysis}} functions.
}

\usage{
  delta_filter(entropy, Smin = 0, Smax = 1)
}

\arguments{
 \item{entropy}{
  An object created by the \code{\link{entropy}} or the \code{\link{dynamic_entropy}}function
 }
 \item{Smin}{
  A value indicating the minimum entropy value. (Smin = 0 by default)
 }
 \item{Smax}{
  A value indicating the maximum entropy value. (Smax = 1 by default)
 }
}

\details{
  The object returned by the \code{\link{entropy}} or the \code{\link{dynamic_entropy}} function contains an entropy score for each element.
  The delta weighting of each element is calculated as follow :
      \deqn{ weighting[i] = \left\{ \begin{array}{ll}  1, & Smin < entropy[i] < Smax \\  0, & otherwise \end{array}  \right. }
}

\value{
  A vector that contains a 0 or 1 weighting score for each element (position in sequence alignment or side chain dihedral angle in trajectory) to limit principal component and top pair analysis to elements within a given entropy range.
}

\author{
  Antoine GARNIER
}

\examples{
  #Importing MSA
  align <- import.msf(system.file("msa/toy_align.msf", package = "Bios2cor"))

  #Creating entropy object
  entropy <- entropy(align)
  
  #Creating delta filter based on entropy
  filter <- delta_filter(entropy, Smin = 0.4, Smax = 0.6)
  
}


