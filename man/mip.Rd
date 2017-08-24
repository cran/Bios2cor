\name{mip}

\Rdversion{1.1}

\alias{mip}

\title{
  MIp(Mutual Information product) function
}

\description{
  Calculates a mutual information score (MI) based on the probability of joint occurrence of events and corrects it with the average product  which is subtracted from the MI score.
}

\usage{
  mip(align, fileHelix= NULL, diag= 0, fileCSV= NULL, gap_val= 0.8, z_score= TRUE)
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
  The MIp score at position [i,j] has been computed with the following formula :
  
  \deqn{{MIp(i,j)} = MI(i,j) - \frac{MI(i,\bar{j})MI(\bar{i},j)}{<MI>}}
  
  with :
    \itemize{
      \item {\emph{\eqn{{MI(i,j) = \sum_{x,y}^{ } p_{x,y}(i,j) ln\frac{p_{x,y}(i,j)}{p_{x}(i)p_{y}(j)}}}}} {}
      \item {\emph{\eqn{MI(i,\bar{j}) = \frac{1}{n-1} \sum_{j \neq i}^{ } MI(i,j)}}} {}
      \item {\emph{\eqn{MI(\bar{i},j) = \frac{1}{n-1} \sum_{i \neq j}^{ } MI(i,j)}}} {}
      \item {\emph{\eqn{<MI> = \frac{2}{n(n-1)} \sum_{i,j}^{ }MI(i,j)}}} {}
    }
    
  and where \emph{\eqn{p_{x,y}(i,j)}} is the frequency of the amino acid pair (x,y) at positions i and j.
      
  N.B. this formula has been widely applied in the field of sequence covariation but favors pairs with high entropy.
}

\value{
  A list of two elements : a matrix containing the correlation scores for each pair of positions and, optionally, a second matrix with the Z-scores 
}

\author{
  Madeline DENIAUD and Antoine GARNIER
}


\references{
 Dunn SD, Wahl LM, Gloor GB. Mutual information without the influence of phylogeny or entropy dramatically improves residue contact prediction. Bioinfor;atics 2008;24:333-340.
 Martin LC, Gloor GB, Dunn SD, Wahl LM. Using infor;ation theory to search for co-evolving residues in proteins. Bioinformatics 2005;21:4116-4124.
}
  
  
\examples{
   align <- import.fasta(system.file("msa/toy2_align.fa", package = "Bios2cor"))

  #Creating MIP object
  mip <- mip(align)
}

\keyword{correlation}
\keyword{MI}
\keyword{mip}
