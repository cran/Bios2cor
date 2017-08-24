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
  entropy(align, fileHelix= NULL)
}

\arguments{
 \item{align}{
  An object created by the \bold{\emph{import.msf}} function which required an alignment file name
 }
 \item{fileHelix}{
  A file that contains the positions of the anchor residues in the sequence alignment. Developped for analysis of GPCR sequence alignments 
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
  Antoine GARNIER
}


\references{
  Shannon CE. A mathematical theory of communication. Bell Syst Techn J 1948;27:379-423.
} 
  
\examples{
  align <- import.msf(system.file("msa/toy_align.msf", package = "Bios2cor"))
 
  #creating entropy object
  #entropy <- entropy(align, fileHelix = helix)
  entropy <- entropy(align)
}

\keyword{entropy}
\keyword{score}
\keyword{Shannon}
