\name{mi}

\Rdversion{1.1}

\alias{mi}

\title{
  Mutual Information (MI) function
}

\description{
  Calculates a mutual information score (MI) based on the probability of joint occurrence of events.
}

\usage{
  mi(align, gap_ratio= 0.2)
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
  The MI score at position [i,j] has been computed with the following formula :
  
    \itemize{
      \item {\emph{\eqn{{MI(i,j) = \sum_{x,y}^{ } p_{x,y}(i,j) ln\frac{p_{x,y}(i,j)}{p_{x}(i)p_{y}(j)}}}}} {}
    }
    
  and where \emph{\eqn{p_{x,y}(i,j)}} is the frequency of the amino acid pair (x,y) at positions i and j.
      
  N.B. this formula has been widely applied in the field of sequence correlation/covariation but favors pairs with high entropy.
}

\value{
  A list of two elements which are numeric matrices containing the MI scores and Z-scores for each pair of elements. 
}

\author{
  Madeline DENIAUD and Marie CHABBERT
}


\references{
 Dunn SD, Wahl LM, Gloor GB. Mutual information without the influence of phylogeny or entropy dramatically improves residue contact prediction. Bioinfor;atics 2008;24:333-340.
 Martin LC, Gloor GB, Dunn SD, Wahl LM. Using infor;ation theory to search for co-evolving residues in proteins. Bioinformatics 2005;21:4116-4124.
}
  
  
\examples{
  #Importing MSA file
  align <- import.fasta(system.file("msa/toy2_align.fa", package = "Bios2cor"))

  #Creating correlation object with MI method for positions with gap ratio < 0.2 (Default)
  mi <- mi(align)
}


