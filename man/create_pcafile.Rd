\name{create_pcafile}

\Rdversion{1.1}

\alias{create_pcafile}

\title{
  Creates a file of coordinates in PCA space 
}

\description{
  Given an object of class 'pca' (result of the \code{\link{centered_pca}} function), stores the coordinates of each element in the PC space in a txt file
}

\usage{
  create_pcafile(
    pca_struct,
    entropy= NULL,
    filepath,
    top_positions= 10,
    top_pca1= NULL,
    top_pca2= NULL
  )
}

\arguments{
 \item{pca_struct}{
  An object created by the \code{\link{centered_pca}} function
 }
 \item{entropy}{
  An object created by the \bold{\emph{entropy}} or the \code{\link{rotamer_entropy}} function 
 }
 \item{filepath}{
  The full path name of the output file where all coordinates on all components are stored
 }
 \item{top_positions}{
  An integer corresponding to the number of elements with highest coordinates on the first or on second principal component that will be save if top_pca1 or top_pca2 are not FALSE. DEFAULT is 10 
 }
 \item{top_pca1}{
  The full path name of the output file where the PC1 and PC2 coordinates of the elements with top PC1 coordinates are saved
 }
 \item{top_pca2}{
  The full path name of the output file where the PC1 and PC2 coordinates of the elements with top PC2 coordinates are saved
 }
}

\details{
  The object returned by the \code{\link{centered_pca}} function contains coordinates in the PC space for each element.
  Each line of the pca file will contain the name of the current element and its coordinates.
  Any line that contains NaN value for X, Y or Z coordinates will be ignored.
  Optionally, the function also returns the top 'top_positions' elements with highest absolute coordinates on PC1 and PC2
}

\value{
  returns a file containing the coordinates of each element in PC space and, optionally, top elements on PC1 or PC2 with their coordinates on PC1 and PC2
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
   
  create_pcafile(pca, entropy = NULL,"pca_results.txt",10,"pc1.txt","pc2.txt")
  
}

\keyword{pca}
\keyword{coordinates}
\keyword{correlation}
