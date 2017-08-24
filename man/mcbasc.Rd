\name{mcbasc}

\Rdversion{1.1}

\alias{mcbasc}

\title{
  McBASC(McLachlan Based Substitution Correlation) function
}

\description{
  Calculates a score for each pair of residus in the sequenec alignment. 
  It relies on a substitution matrix giving a similarity score for each pair of amino acids.
}

\usage{
  mcbasc(align, fileHelix= NULL, diag= 0, fileCSV= NULL, gap_val= 0.8, z_score= TRUE)
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
  A boolean for Z-score normalisation of the covariation matrix. Default is TRUE.
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
  A list of two elements : a matrix containing a correlation score for each pair of rotamer and its normalized version
}

\author{
  Madeline DENIAUD and Antoine GARNIER
}

\examples{
  align <- import.msf(system.file("msa/toy_align.msf", package = "Bios2cor"))

  #Creating McBASC object
  mcbasc <- mcbasc(align)

}

\references{
 (1) Gobel U, Sander C, Schneider R, Valencia A. Correlated mutations and residue contacts in proteins. Proteins 1994;18:309-317.
}

\keyword{correlation}
\keyword{substitution}
\keyword{mcbasc}
