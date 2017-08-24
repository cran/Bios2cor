\name{shuffle_positions}

\Rdversion{1.1}

\alias{shuffle_positions}

\title{
  Amino acids shuffle
}

\description{
 Given an alignment object, created by the \code{\link{import.fasta}} or \code{\link{import.msf}} function, for each position in the alignment, shuffles the residues at this position to create a random alignment, and calculates the correlation matrix with the selected method.  
}

\usage{
  shuffle_positions(
    align,
    method,
    fileHelix= NULL,
    gap_val= 0.8,
    z_score= TRUE,
    nb_iterations= 5
  )
}

\arguments{
 \item{align}{
   An object created by the \bold{\emph{import.msf}} or \bold{\emph{import.fasta}} function which requires an alignment file name
 }
 \item{method}{
   A string corresponding to selected correlation method. This should be one out of "OMES", "MIP","ELSC" and "MCBASC"
 }
 \item{fileHelix}{
   AA file that contains the positions of the anchor residues in the sequence alignment. For use with GPCRs
 }
 \item{gap_val}{
   Authorized gap proportion per position. This value can be between 0 and 0.8 (0.8 by default), which means that positions with more than 80 percent of gaps will not be considered 
 }
 \item{z_score}{
    A logical value to perform a Z-score normalisation of the covariation matrix (TRUE) or not (FALSE). Default is TRUE
 }
 \item{nb_iterations}{
   Number of times that alignment must be shuffled. Default is 5
 }
}

\details{
  At each position, sequences from the initial alignment will be shuffled to yield a random alignment with the same average composition at each position as the initial alignment and the same sequence entropy profile. The correlation matrix of the resulting alignment will be calculated "nb_iterations" times, then averaged. 
}

\value{
   Returns a numeric matrix which is the averaged matrix obtained from re-itering "nb_iterations" times the computation of the correlation matrix from a random alignment with the same composition as the initial alignment
}

\author{
  Antoine GARNIER
}

\examples{
  #Importing sequences alignment
  align <- import.msf(system.file("msa/toy_align.msf", package = "Bios2cor"))

  #Shuffling alignment
  corr <- shuffle_positions(align, method="OMES", fileHelix= NULL, gap_val= 0.8, z_score= TRUE,
       nb_iterations= 3)
}

\keyword{shuffle}
\keyword{align}
\keyword{correlation}
 
