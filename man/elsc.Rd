\name{elsc}

\Rdversion{1.1}

\alias{elsc}

\title{
   Explicit Likelihood of Subset Covariation (ELSC) function
}

\description{
  Calculates a score based on rigorous statistics of correlation/covariation in a perturbation-based algorithm. It measures how many possible subsets of size n would have the composition found in column j in the subset alignment defined by the perturbation in column i, and in the ideal subset
(i.e., in a subset with the amino acid distribution equal to the total alignment).
}

\usage{
  elsc(align, gap_ratio = 0.2)
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
  A list of two elements which are numeric matrices containing the ELSC scores and Z-scores for each pair of elements. 
}

\author{
  Madeline DENIAUD and Marie CHABBERT
}


\references{
 Dekker JP, Fodor A, Aldrich RW, Yellen G. A perturbation-bqsed method for calculating explicit likelihood of evolutionary covariance in multiple sequence alignements. Bioinformatics 2004;20:1565-1572.
 }



\examples{
  #Importing MSA file 
  align <- import.msf(system.file("msa/toy_align.msf", package = "Bios2cor"))

  #Creating correlation object with ELSC method for positions with gap ratio < 0.1
  elsc <- elsc(align, gap_ratio = 0.1)

  #Creating correlation object with ELSC method for positions with gap_ratio < 0.2 (Default)
  #elsc <- elsc(align)
}


