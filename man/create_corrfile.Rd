\name{create_corrfile}

\Rdversion{1.1}

\alias{create_corrfile}

\title{
  Creates a correlation file 
}

\description{
  Given a correlation object, stores the correlation score of each pair of elements in a txt file.
}

\usage{
  create_corrfile(corr_struct, entropy= NULL, auto_pairing= FALSE, filepath)
}

\arguments{
 \item{corr_struct}{
  An object created by a correlation function (\code{\link{omes}}, \code{\link{mip}}, \code{\link{elsc}}, \code{\link{mcbasc}}, \code{\link{rotamer_circular}}, \code{\link{rotamer_omes}}, \code{\link{rotamer_mip}}
 }
 \item{entropy}{
  An object created by the \bold{\emph{entropy}} function 
 }
 \item{auto_pairing}{
   Option that allows/removes correlation scores between dihedral angles of the same residue in the correlation matrix (valid only for side chain correlations in molecular simulations) 
 }
 \item{filepath}{
  The full path name of the output file
 }
}

\details{
  The correlation structure contains a score for each pair of elements [i,j]. Elements represent positions in the sequence alignments and side chain dihedral angles in molecular simulations.
  If entropy = NULL, each line of the correlation file will contain element  i, element  j and score [i,j].
  If entropy is not NULL, each line of the correlation file will contain element  i, element  j, entropy[i], entropy[j] and score [i,j]
}

\value{
  A txt file containing the correlation scores and optionally the entropies
}

\author{
  Antoine GARNIER
}

\examples{
  msf <- system.file("msa/toy_align.msf", package = "Bios2cor")
  align <- import.msf(msf)

  #Creating OMES object
  omes <- omes(align,fileHelix= NULL , diag= 0, fileCSV= NULL, gap_val= 0.8, z_score= TRUE)
  omes <-omes$normalized
  
  create_corrfile(omes, entropy= NULL, auto_pairing= FALSE, "corr_results.txt")
}

\keyword{correlation}
