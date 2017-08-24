\name{create_screeplot}

\Rdversion{1.1}

\alias{create_screeplot}

\title{
  Creates PCA screeplot 
}

\description{
  Given a PCA structure (result of the \code{\link{centered_pca}} function), creates a screeplot of the positive eigen values
}

\usage{
  create_screeplot(pca_struct, pca_index= NULL, filepath)
}

\arguments{
 \item{pca_struct}{
  An object created by the \code{\link{centered_pca}} function
 }
 \item{filepath}{
  The full path name of the output file 
 }
 \item{pca_index}{
  An integer indicating the number of selected eigenvalues. Its value must not exceed correlation matrix length. DEFAULT is NULL (all positive eigen values are plotted)
 }
}

\value{
  A screeplot of selected eigen values 
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

  pca <- centered_pca(omes, m = NULL, pc= NULL, dec_val= 5,eigenvalues_csv= NULL)

  create_screeplot(pca, pca_index= NULL, "~/screeplot.png")
  
}

\keyword{pca}
\keyword{screeplot}
\keyword{variance}
