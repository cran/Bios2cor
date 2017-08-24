\name{corr_contact}

\Rdversion{1.1}

\alias{corr_contact}

\title{
  Calculates the number of contacts for each element of the top pairs 
}

\description{
  Given a correlation object, calculates the number of pairs (contacts) each element in the top X pairs are involved in
}

\usage{
  corr_contact(corr, top= 25, contact_filepath= NULL, score_top_filepath= NULL)
}

\arguments{
 \item{corr}{
  An object created by a correlation function (\code{\link{omes}}, \code{\link{mip}}, \code{\link{elsc}}, \code{\link{mcbasc}}, \code{\link{rotamer_circular}}, \code{\link{rotamer_omes}}, \code{\link{rotamer_mip}}).
 } 
 \item{top}{
  A integer indicating the number of top pairs used for this analysis
 } 
 \item{contact_filepath}{
  The full path name of the output contact file
 }
 \item{score_top_filepath}{
   The full path name of the output top score file
 }
}

\details{
  This function sorts element pairs by correlation score and analyses the top X pairs to determine the number of pairs (contacts) each element of the top X pairs is involved in.
}

\value{
  A vector containing top elements names
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
   
  contact <- corr_contact(omes, top= 25, "contact_file.txt","top_scores.txt")
}

\keyword{contact}
\keyword{top}
\keyword{correlation}
