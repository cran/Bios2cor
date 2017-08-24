\name{import.msf}
\Rdversion{1.1}
\alias{import.msf}
\title{
Reads a multiple sequence alignment file in MSF format
}
\description{
Reads a Multiple Sequence Alignment (MSA) file in MSF format (.msf extension).
}
\usage{
import.msf(file, aa.to.upper = TRUE, gap.to.dash = TRUE, log.file = NULL)
}
\arguments{
  \item{file}{
  a string of characters to indicate the name of the MSA file to be read.
}
  \item{aa.to.upper}{
  a logical value indicating whether amino acids should be converted to upper case (TRUE) or not (FALSE). Default is TRUE.
}
  \item{gap.to.dash}{
  a logical value indicating whether the dot (.) and tilde (\eqn{\sim}) gap symbols should be converted
  to the dash (-) character (TRUE) or not (FALSE). Default is TRUE.
}
  \item{log.file}{
  a file containing errors that occurred when trying to read fasta file
}
}
\details{
Initially, Multiple Sequence Format (MSF) was the multiple sequence alignment format of
the Wisconsin Package (WP) or GCG (Genetic Computer Group). This package is a suite of over 130
sequence analysis programs for database searching, secondary structure prediction
or sequence alignment. Presently, numerous multiple sequence alignment editors
(Jalview and GeneDoc for example) can read and write MSF files.
  
MSF file displays several specificities:
\itemize{
\item a header containing sequence identifiers and characteristics (length, check and weight).
\item a separator symbolized by 2 slashes (//).
\item sequences of identifiers, displayed by consecutive blocks.
}
}
\value{
A object of class 'align', which is a named list whose elements correspond to sequences, in the form of character vectors.
}
\author{
Julien Pele
}
\note{
\code{import.msf} checks the presence of duplicated identifiers in header. Sequences whose
identifiers are missing in header are ignored.
}
\seealso{
\code{read.alignment} function from \code{seqinr} package.\cr
\code{read.GDoc} function from \code{aaMI} package (archived).
}
\examples{
# reading of the multiple sequence alignment of human GPCRs in MSF format:
  aln <- import.msf(system.file("msa/toy_align.msf", package = "Bios2cor"))
}
\keyword{read,import,msf}

