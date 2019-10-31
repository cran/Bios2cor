\name{centered_pca}

\Rdversion{1.1}

\alias{centered_pca}

\title{
  Performs principal component analysis of a correlation/covariation matrix 
}

\description{
  Given a correlation/covariation matrix, performs the principal component analysis of the centered matrix.
}

\usage{
  centered_pca(corr_matrix, filepathroot= NULL, filter = NULL, pc= NULL, dec_val= 5)
}

\arguments{
 \item{corr_matrix}{
  A score or Zscore matrix created by a correlation/covariation function (\code{\link{omes}}, \code{\link{mip}}, \code{\link{elsc}}, \code{\link{mcbasc}}, \code{\link{dynamic_circular}}, \code{\link{dynamic_omes}}, \code{\link{dynamic_mip}})
 }
 \item{filepathroot}{
  The root of the full path name for the csv and png files where eigen values are stored or displayed. Default is NULL (Two "EIGEN.csv" and "EIGEN.png" files are created). If not NULL, the extensions "_EIGEN.csv" and "_EIGEN.png" are added to the filepathroot. 
 }
 \item{filter}{
  A filter calculated by the \code{\link{delta_filter}} function to limit the analysis to elements within a given entropy/dynamic_entropy range. DEFAULT is NULL (no filter is applied and each element has the same weight equal to the inverse of the number of elements).
 }
 \item{pc}{
  A numeric value indicating the number of principal components to be saved. Default is NULL (all the principal components are saved).
 }
 \item{dec_val}{
  A numeric value corresponding to the precision when the \emph{round} function is used. Default is 5.  
 }
}

\details{
  This function performs a principal component analysis of a correlation/covariation matrix after double centering. It is based on the matrix centering algorithm of the \code{mmds.R} function from the \code{Bios2mds} package. The elements have the same weight except when a delta filter is indicated. In this latter case, only the elements allowed by the delta filter are taken into account.
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
  a numeric matrix representing the coordinates of each element of the correlation/covariation matrix in the PCA space
  }
  \item{source}{
  a named list with 2 elements, the correlation/covariation matrix (\code{cor}) and the delta filter vector (\code{filter}) to limit the analysis to elements within a given entropy range.
  }
 returns also two files: a csv file containing eigen values and a png file displaying eigen values.

}



\author{
  Antoine GARNIER and Marie CHABBERT 
}

\examples{
  #File path for output files
  wd <- tempdir()
  #wd <- getwd()
  file <- file.path(wd,"test_seq_pca") 

  #Importing multiple sequence alignment
  msf <- system.file("msa/toy_align.msf", package = "Bios2cor")
  align <- import.msf(msf)

  #Creating correlation object with OMES method
  omes <- omes(align, gap_ratio= 0.2)

  #Creating entropy object
  entropy <- entropy(align, gap_ratio=0.2)

  #Creating delta filter based on entropy
  filter <- delta_filter(entropy, Smin = 0.2, Smax = 0.6)
  
  #Selecting a correlation/covariation matrix
  matrix_omes <-omes$score

  #Creating PCA object for selected matrix and storing eigen values in csv file
  pca <- centered_pca(matrix_omes, filepathroot= file, filter = filter)

}


