\name{write.pca}

\Rdversion{1.1}

\alias{write.pca}

\title{
  Creates a file of coordinates in PCA space 
}

\description{
  Given an object of class 'pca' (result of the \code{\link{centered_pca}} function), stores the coordinates of each element in the PC space in a txt file
}

\usage{
  write.pca(corr_pca,filepathroot=NULL, pc=NULL, entropy= NULL)
}

\arguments{
 \item{corr_pca}{
  An object created by the \code{\link{centered_pca}} function from a correlation/covariation matrix
 }
 \item{filepathroot}{
  The root for the full path name of the output file where all coordinates on all components are stored. Default is NULL (a "PCA_COORD.csv" file is created). If not NULL, the "_PCA_COORD.csv extension is added to the root name.
 }
 \item{pc}{
  An integer corresponding to the number of principal components for which coordinates of the elements are saved. By default, this number corresponds to the number of components with positive eigen values.
 } 
 \item{entropy}{
  An object created by the \code{\link{entropy}} or the \code{\link{dynamic_entropy}} function. DEfault is NULL. 
 }
}

\details{
  The object returned by the \code{\link{centered_pca}} function contains coordinates in the PC space for each element.
  Each line of the pca file will contain the name of the current element and its coordinates.
  Any line that contains Na value for X, Y or Z coordinates will be ignored.
}

\value{
  returns a file containing the coordinates of each element in PC space.
}

\author{
  Antoine GARNIER
}

\examples{
  #File path for output files
  wd <- tempdir()
  #wd <-getwd() 
  file <- file.path(wd,"test_seq4") 

  #Importing MSA file
  msf <- system.file("msa/toy_align.msf", package = "Bios2cor")
  align <- import.msf(msf)

  #Creating OMES correlation object and selecting correlation matrix
  omes <- omes(align, gap_ratio = 0.2)
  omes <-omes$Zscore

  #Creating PCA object for selected matrix
  pca <- centered_pca(omes, filepathroot= file, filter = NULL, pc= NULL, dec_val= 5)
   
  #Saving coordinates of elements in csv file 
  write.pca(pca, file, pc = 10, entropy = NULL)
  
}


