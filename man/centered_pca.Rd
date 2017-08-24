\name{centered_pca}

\Rdversion{1.1}

\alias{centered_pca}

\title{
  Principal component analysis of a correlation/covariation matrix
}

\description{
  Given a correlation matrix, performs a principal component analysis of a correlation/covariation matrix.
}

\usage{
  centered_pca(mat, m = NULL, pc= 3, dec_val= 3, eigenvalues_csv= NULL)
}

\arguments{
 \item{mat}{
  A matrix created by a correlation function (\code{\link{omes}}, \code{\link{mip}}, \code{\link{elsc}}, \code{\link{mcbasc}}, \code{\link{rotamer_circular}}, \code{\link{rotamer_omes}}, \code{\link{rotamer_mip}})
 }
 \item{m}{
  A weighting filter calculated by a weighting function (\code{\link{gauss_weighting}}, \code{\link{sigmoid_weighting}} or \code{\link{delta_weighting}}). DEFAULT is NULL (no filter is applied and each element has the same weight equal to the inverse of the number of elements).
 }
 \item{pc}{
  A numeric value indicating the number of principal components to be saved.
 }
 \item{dec_val}{
  A numeric value corresponding to the precision when the \emph{round} function is used.
 }
 \item{eigenvalues_csv}{
  A full path name for the csv file where eigen values are stored. Default is NULL (csv file is not created). If not NULL, a png file will also be created with the .csv extension replaced by .png
 }
}

\details{
  This function performs a principal component analysis of a correlation/covariation matrix after double centering. It is based on the matrix centering of the \emph{mmds.R} function from the \emph{Bios2mds} package.The elements may have the same weight or different weights. In this case, a weighting filter, corresponding to the relative weight of each matrix element, is passed as a parameter The weighting filter is calculated by a weighting function based on the element entropy.
}

\value{
  returns an object of class 'pca' which is a named list of four elements: 
  \item{eigen}{
  a numeric vector of the eigenvalues
  }
  \item{eigen.perc}{
  a numeric vector of the relative eigenvalues (eigenvalues divided by the sum of the absolute eigenvalues)
  }
  \item{coord}{
  a numeric matrix representing the coordinates of each element of the correlation matrix in the PCA space
  }
  \item{source}{
  a named list with 2 elements, the correlation matrix (\bold{cor}) and the vector of mass \bold{m} which is used to give a different weight to each element.
  }
}



\author{
  Antoine GARNIER 
}

\examples{
   msf <- system.file("msa/toy_align.msf", package = "Bios2cor")
  align <- import.msf(msf)

  #Creating OMES object
  omes <- omes(align,fileHelix= NULL , diag= 0, fileCSV= NULL, gap_val= 0.8, z_score= TRUE)

  #Creating ENTROPY object
  entropy <- entropy(align)

  #Creating weighting filter
  filter <- gauss_weighting(entropy, L= 0.1)
  
  # Creating PCA structures for OMES method and storing in txt file
  omes <-omes$normalized
  pca <- centered_pca(omes, m= filter, pc= NULL, dec_val= 5,eigenvalues_csv= NULL)
  #create_pcafile(pca, "pca_results.txt",10,"pc1.txt","pc2.txt")
}

\keyword{pca}
\keyword{weighting}
\keyword{correlation}
