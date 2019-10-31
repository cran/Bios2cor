\name{random.msa}
\Rdversion{1.1}
\alias{random.msa}
\title{
Random Alignment
}
\description{
Builds a multiple sequence alignment (MSA) of random sequences.
}
\usage{
random.msa(nb.seq = 100, id = "SEQ", nb.pos = 100, gap = FALSE,
aa.strict = FALSE, align = NULL, align.replace = TRUE)
}
\arguments{
  \item{nb.seq}{
  a numeric value indicating the number of sequences in the random MSA. Default is 100. 
}
\item{id}{
  a string of characters used to tag each sequence name. Default is "SEQ". An
  incremented number is attached to this tag to name each sequence.
}
  \item{nb.pos}{
  a numeric value indicating the length of each sequence in the random MSA. Default is 100.
}
  \item{gap}{
  a logical value indicating whether the gap character should be considered as a supplementary symbol
  (TRUE) or not (FALSE). Default is FALSE. 
}
  \item{aa.strict}{
  a logical value indicating whether only strict amino acids should be taken into account (TRUE)
  or not (FALSE). Default is FALSE.
}
  \item{align}{
  an object of class 'align', obtained from \code{\link{import.fasta}} or \code{\link{import.msf}} function.
  If this parameter is not NULL, the composition of the output sequences is based on the composition of the input sequences. Default is NULL.
}
  \item{align.replace}{
  a logical value indicating random drawing with replacement (TRUE) or without replacement (FALSE) of characters present in \code{align}. Default is FALSE.
}
}
\details{
\code{random.msa} may be used to compare a reference MSA to a random MSA. The random MSA must have the same characteristics as the reference MSA (same number of sequences of same length). 

A procedure can be applied to the random MSA to assess the amount of variance due to random mutations in the reference MSA. 

The \code{subset} function is used for random selection of the amino acids. If a truly random procedure is needed, see \code{random} package.
}
\value{
A named list whose objects correspond to random sequences.
}
\note{
This function has been initially developped in the \code{bios2mds} R package (Julien PELE [aut], Marie CHABBERT [cre]).
}
\author{
Julien PELE
}
\references{
For an application of random MSA see : 

Pele J, Abdi H, Moreau M, Thybert D and Chabbert M (2011) Multidimensional scaling reveals the main evolutionary pathways of class A G-protein-coupled receptors. \emph{PLoS ONE} \bold{6}: e19094. doi:10.1371.
}
\seealso{
\code{permutation} and \code{synsequence} functions from \code{seqinr} package.
}
\examples{
  #Importing MSA file
  aln <- import.fasta(system.file("msa/toy2_align.fa", package = "Bios2cor"))

  #Generating a random sequence alignment with the same characterics as input file
  nb.seq <- length(aln)
  nb.pos <- length(aln[[1]])
  aln.random <- random.msa(nb.seq = nb.seq, nb.pos = nb.pos)
}


