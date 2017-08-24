\name{pca_2d}

\Rdversion{1.1}

\alias{pca_2d}

\title{
  PCA projection on two dimensions
}

\description{
  Given an object of class 'pca'  (result of the \code{\link{centered_pca}} function), draws a graph of the projection of the elements on two dimensions (first and second components by default)
}

\usage{
  pca_2d(pca_struct, abs= 1, ord= 2, filepath)
}

\arguments{
 \item{pca_struct}{
  An object created by the \code{\link{centered_pca}} function
 }
 \item{abs}{
  An integer which corresponds to the component that will be projected onto the x axis. Default is 1 (first component)
 }
 \item{ord}{
  An integer which corresponds to the component that will be projected onto the y axis. Default is 2 (second component)
 }
 \item{filepath}{
  The full path name of the png file that will be creared
 }
}

\value{
  A 2D graph of the projection of the elements on selected two principal components
}

\author{
  Antoine GARNIER
}

\examples{
  msf <- system.file("msa/toy_align.msf", package = "Bios2cor")
  align <- import.msf(msf)

  #Creating OMES object
  omes <- omes(align)

  #Creating ENTROPY object
  entropy <- entropy(align)

  # Creating PCA structures for OMES method and storing in txt file
  omes <-omes$normalized
  pca <- centered_pca(omes, m= NULL, pc= NULL, dec_val= 5,eigenvalues_csv= NULL)
  pca_2d(pca, abs= 1, ord= 2, "pca_2d.png")
  
}

\keyword{pca}
\keyword{2d}
