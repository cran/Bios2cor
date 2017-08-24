\name{create_network}

\Rdversion{1.1}

\alias{create_network}

\title{
  Creates network structure of top elements
}

\description{
  Given a top_contact structure (result of the \code{\link{corr_contact}} function), creates a network to visualize the elements involved in the top scoring pairs and their links
}

\usage{
  create_network(top_contact_struct, filepath)
}

\arguments{
 \item{top_contact_struct}{
  An object created by the \code{\link{corr_contact}} function
 }
 \item{filepath}{
  The full path name of the output file
 }
}

\value{
  A network representing links between elements in the top scoring pairs
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

  create_network(contact, "network_file.pdf")
  
}

\keyword{contact}
\keyword{network}
