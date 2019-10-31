\name{network.plot}

\Rdversion{1.1}

\alias{network.plot}

\title{
  Creates network structure of top elements
}

\description{
  Given a top_pairs object (result of the \code{\link{top_pairs_analysis}} function), creates a network to visualize the elements involved in the top scoring pairs and their links
}

\usage{
  network.plot(top_pairs, filepathroot=NULL)
}

\arguments{
 \item{top_pairs}{
  An object of class 'top_pairs' created by the \code{\link{top_pairs_analysis}} function
 }
 \item{filepathroot}{
  The root of the full path name for the output file. Default is NULL (a "NETWORK.pdf" file will be created). If not NULL, the filepathroot will have the "_NETWORK.pdf" extension. 
 }
}

\value{
  A network representing links between elements in the top scoring pairs
}

\author{
  Antoine GARNIER
}

\examples{
  #File path for output file
  wd <- tempdir()
  #wd <-getwd() 
  file <- file.path(wd,"test_omes") 

  #Importing MSA file
  msf <- system.file("msa/toy_align.msf", package = "Bios2cor")
  align <- import.msf(msf)

  #Creating OMES correlation object for positions with gap ratio < 0.2
  omes <- omes(align, gap_ratio= 0.2)
  
  #Selecting a correlation matrix
  omes <-omes$Zscore
   
  #Analyzing pairs with top scores and creating 'top_pairs' object
  top_pairs <- top_pairs_analysis(omes, top = 25, file)

  #Plotting the network structure of top pairs in pdf file 
  network.plot(top_pairs, file)
}


