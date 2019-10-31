\name{write.entropy}

\Rdversion{1.1}

\alias{write.entropy}

\title{
  Writes and displays entropy values
}

\description{
  Given an entropy/dynamic_entropy object, writes each element (position or dihedral angle) and its entropy/dynamic_entropy value in a csv file, and displays the histogram. 
}

\usage{
  write.entropy(entropy, filepathroot = NULL)
}

\arguments{
 \item{entropy}{
  An object created by the \code{\link{entropy}} or the \code{\link{dynamic_entropy}} function.
 }
 \item{filepathroot}{
  The root of the full path name for the entropy output file. Default is NULL (a "ENTROPY.csv" file is created). If filepathroot is not NULL, a "_ENTROPY.csv" extension is added to the file path root. 
 }
}

\details{
  The \code{\link{entropy}} function calculates entropy score for each position of an alignment file. The \code{\link{dynamic_entropy}}function calculate a "dynamic entropy" score for each side chain dihedral angle of a protein during a molecular simulations. 
}


\value{
  A csv file containing the elements and their scores. A png file displaying the histogram of the scores.
}

\author{
  Antoine GARNIER
}

\examples{
  #File path for output files
  wd <- tempdir()
  #wd <-getwd() 
  file <- file.path(wd,"test_seq") 

  #Importing multiple sequence alignment
  align <- import.msf(system.file("msa/human_gpcr.msf", package = "Bios2cor"))

  #Creating entropy object
  entropy <- entropy(align, gap_ratio = 0.2)
  
  #Writing results to csv file
  write.entropy(entropy, file)
}


 
