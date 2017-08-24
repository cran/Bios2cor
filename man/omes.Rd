\name{omes}

\Rdversion{1.1}

\alias{omes}

\title{
  OMES(Observed minus Expected Squared) function
}

\description{
  Calculates the difference between the observed and expected occurrences of each possible pair of amino acids (x, y) at positions i and j of the alignment
}

\usage{
  omes(align, fileHelix= NULL , diag= 0, fileCSV= NULL, gap_val= 0.8, z_score= TRUE)
}

\arguments{
 \item{align}{
  An object of class 'align' created by the \bold{\emph{import.msf}} or the \bold{\emph{import.fasta}} function from a sequence alignment
 }
 \item{fileHelix}{
   A string of characters that indicates the file containing the positions of the anchor residues in the sequence alignment. To be used for the analysis of GPCR sequences. Default is NULL.
 }
 \item{diag}{
  A numeric value indicating the score of the diagonal elements in the scoring matrix. Default is 0.
 }
 \item{fileCSV}{
  A string of characters indicating the name of the csv file where the output matrix will be saved. Default is NULL.
 }
 \item{gap_val}{
   Numeric value indicating the gap ratio at a given position for this position to be taken into account. This value must be between 0 and 0.8. Default is 0.8, which means that positions with more than 80 percent of gaps will not be taken into account.
 }
 \item{z_score}{
  A logical value to perform a Z-score normalisation of the covariation matrix (TRUE) or not (FALSE). Default is TRUE.
 }
}

\details{
  The OMES score at position [i,j] has been computed with the following formula :
  
  \deqn{{OMES(i,j)} = \frac{1}{N(i,j)} \sum_{x,y}^{ }(N_{x,y}^{obs}(i,j)-N_{x,y}^{ex}(i,j))^2}
  
  with : \emph{\eqn{N_{x,y}^{ex}(i,j) = p_{x}(i)p_{y}(j)N}}

  where :
    \itemize{
      \item {\emph{\eqn{N_{x,y}^{obs}(i,j)}}}{is number of times that each (x,y) pair is observed at positions i and j}
      \item {\emph{\eqn{N_{x,y}^{ex}(i,j)}}}{is number of times that each (x,y) pair would be expected, based on the frequency of residues x and y at positions i and j}
      \item {\emph{\eqn{N(i,j)}}}{is the number of sequences in the alignment with non-gapped residues at positions i and j}
      \item {\emph{\eqn{p_{x}(i)}}}{is the frequency of amino acid x at position i}
      \item {\emph{\eqn{p_{y}(j)}}}{is the frequency of amino acid y at position j}
    }
}

\value{
  A list of two elements : a matrix containing a correlation score for each pair of elements and, optionally, a matrix containing the Z-scores
}

\author{
  Madeline DENIAUD and Antoine GARNIER
}

\references{
Fodor AA and Aldrich RW. Influence of conservation on calculations of amino acid covariance in multiple
sequence alignments. Proteins. 2004;56(2):211-21.
}

\examples{
   msf <- system.file("msa/toy_align.msf", package = "Bios2cor")
  align <- import.msf(msf)

  #Creating OMES object
  #omes <- omes(align,fileHelix= NULL , diag= 0, fileCSV= NULL, gap_val= 0.8, z_score= TRUE)
  omes <- omes(align)
}

\keyword{correlation}
\keyword{homogeneity}
\keyword{omes}
