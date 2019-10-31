\name{pca_screeplot}

\Rdversion{1.1}

\alias{pca_screeplot}

\title{
  Creates PCA screeplot 
}

\description{
  Given a PCA structure (result of the \code{\link{centered_pca}} function), creates a screeplot of the positive eigen values
}

\usage{
  pca_screeplot(pca_struct, filepathroot=NULL)
}

\arguments{
 \item{pca_struct}{
  An object created by the \code{\link{centered_pca}} function
 }
 \item{filepathroot}{
  The root of the full path name for the output file. Default is NULL (a "EIGEN.png" file will be created). If not NULL, the filepathroot will have the "_EIGEN.png" extension. 
 }
}

\value{
  A screeplot of positive eigen values 
}

\author{
  Antoine GARNIER
}

\examples{
  #File path for output files
  wd <- tempdir()
  #wd <-getwd() 
  file <- file.path(wd,"test_seq3") 

  #Importing MSA file
  msf <- system.file("msa/toy_align.msf", package = "Bios2cor")
  align <- import.msf(msf)

  #Creating OMES correlation object
  omes <- omes(align, gap_ratio = 0.2)

  #Selecting correlation matrix
  omes <-omes$Zscore

  #Creating PCA object for selected correlation matrix and saving eigen values
  pca <- centered_pca(omes, filepathroot= file, pc= NULL, dec_val= 5, filter = NULL)

  #Plotting scree plot 
  pca_screeplot(pca, file)
  #pca_screeplot(pca)  
}


