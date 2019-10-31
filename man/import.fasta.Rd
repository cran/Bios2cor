\name{import.fasta}
\Rdversion{1.1}
\alias{import.fasta}
\title{
Reads a file in FASTA format
}
\description{
 Reads a Multiple Sequence Alignment (MSA) file in FASTA format (.fasta or .fa extension).
}
\usage{
 import.fasta(file, aa.to.upper = TRUE, gap.to.dash = TRUE, log.file = NULL)
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
    to dash (-) character (TRUE) or not (FALSE). Default is TRUE.
  }
  \item{log.file}{
    a file containing errors that occurred when trying to read fasta file
  }
}

\details{
  Initially, FASTA (for FAST-ALL) was the input format of the FASTA program, used for protein comparison and searching in databases.
  Presently, FASTA format is a standard format for biological sequences.
  
  The FASTA formatted file of a single sequence displays:
  \itemize{
    \item a single-line description beginning with a greater-than (>) symbol. The following word is the identifier.
    \item followed by any number of lines, representing biological sequence.
  }
  For multiple alignments, the FASTA formatted sequences are concatenated to create a multiple
  FASTA format.
}

\value{
  A object of class 'align', which is a named list whose elements correspond to sequences, in the form of character vectors.
}

\references{
  Pearson WR and Lipman DJ (1988) Improved tools for biological sequence comparison.
  \emph{Proc Natl Acad Sci U S A} \bold{27}:2444-2448.
}

\note{
  The \code{import.fasta} function was developed for the \code{bios2mds} R package (Julien PELE [aut], Jean-Michel BECU [aut], Marie CHABBERT [cre]).
.
}

\author{
  Julien PELE
}


\examples{
  # Importing MSA file in FASTA format
  aln <- import.fasta(system.file("msa/toy2_align.fa", package = "Bios2cor"))
}



