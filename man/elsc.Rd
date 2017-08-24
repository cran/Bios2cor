\name{elsc}

\Rdversion{1.1}

\alias{elsc}

\title{
  ELSC(Explicit Likelihood of Subset Covariation) function
}

\description{
  Calculates a score based on rigorous statistics of covariation in a perturbation-based algorithm. It measures how many possible subsets of size n would have the composition found in column j in the subset alignment defined by the perturbation in column i, and in the ideal subset
(i.e., in a subset with the amino acid distribution equal to the total alignment).
}

\usage{
  elsc(
    align,
    fileHelix= NULL,
    diag= 0,
    fileCSV= NULL,
    gap_val= 0.8,
    double_passing= FALSE,
    z_score= TRUE
  )
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
 \item{double_passing}{
   Boolean to calculate correlation score twice : once from first position to last position then from last to first. Results are summed then divided by 2. DEfault is FALSE.
 }
 \item{z_score}{
   A boolean for Z-score normalisation of the covariation matrix. Default is TRUE.
 }
}

\details{
  The ELSC score at position [i,j] has been computed with the following formula :
  
  \deqn{ELSC(i,j) = -ln\prod_{y}^{ } \frac{{{N_{y(j)}}\choose{n_{y(j)}}}}{{{N_{y(j)}}\choose{m_{y(j)}}}}}
  
  As a reminder, a binomial coefficient \emph{\eqn{{N}\choose{k}}} is computed as follow :
  \emph{\deqn{{{N}\choose{k}} = \frac{N!}{k!(N-k)!}}}
  
  where :
    \itemize{
      \item {\emph{\eqn{N_{y(j)}}}} {is the number of residues y at position j in the total (unperturbed) sequence alignment}
      \item {\emph{\eqn{n_{y(j)}}}} {is the number of residues y at position j in the subset alignment defined by the perturbation in column i}
      \item {\emph{\eqn{m_{y(j)}}}} {is the number of residues y at position j in the ideal subset (i.e., in a subset with the amino acid distribution equal to the total alignment)}
    }
}

\value{
  A list of two elements : a matrix containing the ELSC scores for each pair of elements and, optionally, a matrix containing the Z-scores
}

\author{
  Madeline DENIAUD and Antoine GARNIER
}


\references{
 Dekker JP, Fodor A, Aldrich RW, Yellen G. A perturbation-bqsed method for calculating explicit likelihood of evolutionary covariance in multiple sequence alignements. Bioinformatics 2004;20:1565-1572.
 }



\examples{
  align <- import.msf(system.file("msa/toy_align.msf", package = "Bios2cor"))

  #Creating ELSC object
  elsc <- elsc(align)

}

\keyword{correlation}
\keyword{perturbation}
\keyword{deviation}
\keyword{elsc}
