\name{mcbasc}

\Rdversion{1.1}

\alias{mcbasc}

\title{
  McBASC (McLachlan Based Substitution Correlation) function
}

\description{
  Calculates a score for each pair of residus in the sequenec alignment. 
  It relies on a substitution matrix giving a similarity score for each pair of amino acids.
}

\usage{
  mcbasc(align, gap_ratio= 0.2)
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
  The McBASC score at position [i,j] has been computed with a formula which was initially proposed by Valencia and coworkers(1) as follow :
  
  \deqn{{McBASC(i,j)} = \frac{1}{N^2\sigma(i)\sigma(j)} \sum_{k,l}^{ } (SC_{k,l}(i)-SC(i))(SC_{k,l}(j)-SC(j))}
  
  where :
    \itemize{
      \item {\emph{\eqn{SC_{k,l}(i)}}} {is the score for the amino acid pair present in sequences k and l at position i}
      \item {\emph{\eqn{SC_{k,l}(j)}}} {is the score for the amino acid pair present in sequences k and l at position j}
      \item {\emph{\eqn{SC(i)}}} {is the average of all the scores \emph{\eqn{SC_{k,l}(i)}}}
      \item {\emph{\eqn{SC(j)}}} {is the average of all the scores \emph{\eqn{SC_{k,l}(j)}}}
      \item {\emph{\eqn{\sigma(i)}}} {is the standard deviation of all the scores \emph{\eqn{SC_{k,l}(i)}}}
      \item {\emph{\eqn{\sigma(j)}}} {is the standard deviation of all the scores \emph{\eqn{SC_{k,l}(j)}}}
    }
}

\value{
  A list of two elements which are numeric matrices containing the mcbasc scores and Z-scores for each pair of elements.
}

\author{
  Madeline DENIAUD and Marie CHABBERT
}

\references{
 (1) Gobel U, Sander C, Schneider R, Valencia A. Correlated mutations and residue contacts in proteins. Proteins 1994;18:309-317.
}


\examples{
  #Importing MSA file
  align <- import.msf(system.file("msa/toy_align.msf", package = "Bios2cor"))

  #Creating correlation object with McBASC method for positions with gap_ratio < 0.2 (Default)
  mcbasc <- mcbasc(align, gap_ratio = 0.2)
}



