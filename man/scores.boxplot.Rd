\name{scores.boxplot}

\Rdversion{1.1}

\alias{scores.boxplot}

\title{
  Creates boxplots of correlation/covariation scores 
}

\description{
  Given a list of correlation/covariation matrices, build boxplots for comparative purposes.
}

\usage{
  scores.boxplot(corr_matrix_list, name_list, filepathroot=NULL, elite=25, high=275)
}

\arguments{
 \item{corr_matrix_list}{
  A list of correlation/covariation matrices to be compared
 }
 \item{name_list}{
  The names of the correlation/covariation matrices 
 }
 \item{filepathroot}{
  The root of the full path name for the output file. Default is NULL (a BOXPLOT.png file will be created). If not NULL, the "_BOXPLOT.png" extension is added to the filepathroot. 
 }
 \item{elite}{
  An integer to determine the number of pairs with the highest and lowest scores (e.g. 25: pairs ranked 1 to 25 in decreasing or increasing order) to be colored with the "elite" color codes. Default is 25.
 }
 \item{high}{
  An integer to determine the number of pairs with the next highest and lowest scores (e.g. 275: pairs ranked 26 to 275 in decreasing or increasing order) to be colored with the "high" color codes. Default is 275.
 }
}

\details{
  The correlation/covariation matrices contain the correlation/covariation scores for each pair of elements [i,j].
  The boxplots will allow comparing these scores using color codes : the highest values are dark blue, the next highest values are light blue, the lowest values are red and the next lowest values are pink.
}

\value{
  A pdf figure with boxplots to compare correlation/covariation scores
}

\references{
For an application of these boxplots, see : 

Pele J, Abdi H, Moreau M, Thybert D and Chabbert M (2011) Multidimensional scaling reveals the main evolutionary pathways of class A G-protein-coupled receptors. \emph{PLoS ONE} \bold{6}: e19094. doi:10.1371.
}


\author{
  Julien PELE and Antoine GARNIER
}

\examples{
  #File path for output file
  wd <- tempdir()
  #wd <-getwd() 
  file <- file.path(wd,"test_seq") 

  #Importing MSA file
  msf <- system.file("msa/toy_align.msf", package = "Bios2cor")
  align <- import.msf(msf)

  #Creating OMES correlation object
  omes <- omes(align, gap_ratio = 0.2)

  #Creating MIP correlation object
  #mip <- mip(align, gap_ratio = 0.2)

  #Selecting correlation matrices
  omes <-omes$Zscore
  #mip <- mip$Zscore

  #Creating a list of matrices and plotting the boxplots in a graph
  #Two matrices
  #corr_matrix_list <- list(omes,mip)
  #name <- c("omes","mip")
  #One matrix
  corr_matrix_list <- list(omes)
  name <- c("omes")
  scores.boxplot(corr_matrix_list, name, file, 25, 275)
}


