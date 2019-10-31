\name{entropy}

\Rdversion{1.1}

\alias{entropy}

\title{
Entropy score
}

\description{
  Measures the entropy score of each position in a sequence alignment
}

\usage{
  entropy(align, gap_ratio=0.2)
}

\arguments{
 \item{align}{
  An object created by the \code{\link{import.msf}} or the \code{\link{import.fasta}} function from a multiple sequence alignment file
 }
 \item{gap_ratio}{
  Numeric value between 0 and 1 indicating the maximal gap ratio at a given position in the MSA for this position to be taken into account. 1 is excluded (positions with gaps only). Default is 0.2, positions with more than 20 percent of gaps will not be taken into account in the analysis. When gap_ratio is 1 or close to 1, only positions with at least 1 aa are taken into account (positions with only gaps are excluded). 
 }

}

\details{
  The entropy score S at position i has been computed with a formula derived from the Shannon's entropy as follow : 
  \deqn{{S(i)} = - \sum_{x}^{ } {p} _{x}(i) \log_{20} {p} _{x}(i)}

  where :
    \itemize{
      \item {\emph{i}} {is the position in the sequence}
      \item {\emph{x}} {is the sequence index}
      \item {\emph{ \eqn{{p} _{x}(i)}  }} {represents the frequency of residue x at position i}
    }
}

\value{
  A vector containing an entropy value for each position in the alignment
}

\author{
  Antoine GARNIER and Marie CHABBERT
}


\references{
  Shannon CE. A mathematical theory of communication. Bell Syst Techn J 1948;27:379-423.
} 
  
\examples{
  #Importing MSA file
  align <- import.msf(system.file("msa/human_gpcr.msf", package = "Bios2cor"))

  #creating entropy object for positions with gap ratio < 0.5
  entropy <- entropy(align,gap_ratio=0.5)
}


