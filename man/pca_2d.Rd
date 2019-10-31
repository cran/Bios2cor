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
  pca_2d(pca_struct, abs= 1, ord= 2, filepathroot=NULL)
}

\arguments{
 \item{pca_struct}{
  An object created by the \code{\link{centered_pca}} function
 }
 \item{abs}{
  An integer which corresponds to the x axis of the projection plane. Default is 1 (first component)
 }
 \item{ord}{
  An integer which corresponds to the y axis of the projection plane. Default is 2 (second component)
 }
 \item{filepathroot}{
  The root of the full path name for the output file. Default is NULL (a "PCA_abs_ord.png" file will be created). If not NULL, the filepathroot will have the "_PCA_abs_ord.png" extension. 
 }
}

\value{
  A 2D graph of the projection of the elements on selected two principal components
}

\author{
  Antoine GARNIER
}

\examples{
  #File path for output files
  wd <- tempdir()
  #wd <-getwd() 
  file <- file.path(wd,"test_seq2") 

  #Importing MSA file
  msf <- system.file("msa/toy_align.msf", package = "Bios2cor")
  align <- import.msf(msf)

  #Creating OMES correlation object
  omes <- omes(align, gap_ratio = 0.2)

  #Creating entropy object
  #entropy <- entropy(align)

  #Selecting correlation matrix
  omes <-omes$Zscore

  #Creating PCA object for selected correlation matrix and saving eigen values
  pca <- centered_pca(omes, filepathroot= file, pc= NULL, dec_val= 5, filter= NULL)

  #Creating 2D plot of elements projected on selected [abs,ord] plane 
  pca_2d(pca, abs = 1, ord = 3, file)

  #Creating 2D plot of elements projected on default [1,2] plane 
  #pca_2d(pca)  
}


