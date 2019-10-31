\name{omes}

\Rdversion{1.1}

\alias{omes}

\title{
  OMES (Observed minus Expected Squared) function
}

\description{
  Calculates the difference between the observed and expected occurrences of each possible pair of amino acids (x, y) at positions i and j of the alignment
}

\usage{
  omes(align, gap_ratio = 0.2)
}

\arguments{
 \item{align}{
  An object of class 'align' created by the \code{\link{import.msf}} or the \code{\link{import.fasta}} function from a sequence alignment
 }
 \item{gap_ratio}{
  Numeric value between 0 and 1 indicating the maximal gap ratio at a given position in the MSA for this position to be taken into account. Default is 0.2, positions with more than 20 percent of gaps will not be taken into account in the analysis. When gap_ratio is 1 or close to 1, only positions with at least 1 aa are taken into account (positions with only gaps are excluded). 
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
  A list of two elements which are numerical matrices containing the OMES scores and Z-scores for each pair of elements. 
}

\author{
  Jean-Miche BECU, Madeline DENIAUD, and Marie CHABBERT
}

\references{
Fodor AA and Aldrich RW. Influence of conservation on calculations of amino acid covariance in multiple
sequence alignments. Proteins. 2004;56(2):211-21.
}

\examples{
  #Importing MSA file
  msf <- system.file("msa/toy_align.msf", package = "Bios2cor")
  align <- import.msf(msf)

  #Creating correlation object with OMES method for positions with gap ratio < 0.2 (Default)
  omes <- omes(align, gap_ratio = 0.2)
  #omes <- omes(align)
}


