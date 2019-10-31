\name{top_pairs_analysis}

\Rdversion{1.1}

\alias{top_pairs_analysis}

\title{
  Calculates the number of contacts for each element of the top pairs 
}

\description{
  Given a correlation object, calculates the number of pairs (contacts) each element in the top X pairs are involved in
}

\usage{
  top_pairs_analysis(corr_matrix, filepathroot=NULL, top=25, entropy=NULL, filter=NULL)
}

\arguments{
 \item{corr_matrix}{
  One of the matrices created by a correlation/covariation function (\code{\link{omes}}, \code{\link{mip}}, \code{\link{elsc}}, \code{\link{mcbasc}}, \code{\link{dynamic_circular}}, \code{\link{dynamic_omes}}, \code{\link{dynamic_mip}}).
 } 
 \item{filepathroot}{
  The root of the full path names of the output files for top_pairs_analysis. Default is NULL (two csv files are created : TOPn_CONTACTS.csv and TOPn_SCORES.csv, where n is the number of top pairs). If not NULL, extentions "-TOPn_CONTACTS.csv" and "_TOPn_SCORES.csv" are added to the root name.
 }
 \item{top}{
  A integer indicating the number of top pairs used for this analysis. Default is 25.
 } 
 \item{entropy}{
  An object created by the entropy or dynamic_entropy function. Default is NULL.
 } 
 \item{filter}{
   A vector created by the \code{\link{delta_filter}} function to limit the analysis to elements in the given entropy/dynamic_entropy range. DEfault is NULL.
 } 
}

\details{
  This function sorts element pairs by correlation/covariation scores and analyzes the top X pairs to determine the number of pairs (contacts) each element of the top X pairs is involved in. If filter is TRUE, only the scores of elements in the delta filter defined entropy range are taken into account. Results are written as .csv files.
}

\value{
returns an object of class 'top_pairs' which is a named list of four elements for subsequent network representation with the \code{network.plot} function:
  \item{pair_i}{
  a vector containing the name of element i in the ordered top pairs
  }
  \item{pair_j}{
  a vector containing the name of element j in the ordered top pairs
  }
  \item{positions}{
  a vector containing the positions in the top n pairs
  }
  \item{contacts}{
  a vector containing the number of contacts of the positions in the top n pairs
  }
returns also two .csv files containing scores and contacts of the top n pairs for subsequent network representation with Cytoscape.

}

\author{
  Antoine GARNIER and Marie CHABBERT
}


\examples{
  #File path for output files
  wd <- tempdir()
  #wd <-getwd() 
  file <- file.path(wd,"test_seq") 

  #Importing MSA file
  msf <- system.file("msa/toy_align.msf", package = "Bios2cor")
  align <- import.msf(msf)

  #Creating entropy object
  entropy <- entropy(align)

  #Creating OMES correlation object
  omes <- omes(align, gap_ratio = 0.2)

  #Selecting correlation matrix
  omes <-omes$Zscore
   
  #Creating top_pairs object and writing scores and contacts to csv files  
  top_pairs <- top_pairs_analysis(omes, file, top = 25, entropy=entropy)
}


