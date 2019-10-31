\name{scores_entropy.plot}

\Rdversion{1.1}

\alias{scores_entropy.plot}

\title{
 Scores versus entropy graph
}

\description{
  The aim of this graph is the visualization of the entropy values of the elements in top and bottom pairs. 
}

\usage{
  scores_entropy.plot(entropy, corr_matrix, 
    filepathroot=NULL, 
    elite=25, 
    high=275, 
    filter=NULL)
}

\arguments{
 \item{entropy}{
  An object created by the \code{\link{entropy}} (sequence alignment) or the \code{\link{dynamic_entropy}} (trajectory) function. 
 }
 \item{corr_matrix}{
  A correlation/covariation matrix created by one of the correlation/covariation functions: 
 }
  \item{filepathroot}{
  The root of the full path name of the graph that will be created. Default is NULL (a "ei_ej.pdf" file will be created). If not NULL, the filepathroot will have the "_ei_ej.pdf" extension. 
 }
 \item{elite}{
  An integer to determine the number of pairs with the highest and lowest scores (e.g. 25: pairs ranked 1 to 25 in decreasing or increasing order) to be colored with the "elite" color codes. Default is 25.
 }
 \item{high}{
  An integer to determine the number of pairs with the next highest and lowest scores (e.g. 275: pairs ranked 26 to 275 in decreasing or increasing order) to be colored with the "high" color codes. Default is 275.
 }
 \item{filter}{
  A vector created by the \code{\link{delta_filter}} function, to limit the analysis to elements within a given entropy range. When filter is not NULL, only top elite and high pairs are visualized. Default is NULL. 
 }
}

\details{
  Using the result of a correlation/covariation method and an entropy structure, creates a graph comparing correlation/covariation scores with entropy values. Each pair of elements (i,j) is placed in the graph with (entropy[i] ; entropy[j]) as coordinates. In the absence of filter, the color code of each point is based on its correlation/covariation score (dark and light blue for top elite and high values, red and pink for bottom elite and high values). In the presence of an entropy based filter, only top elite and high scores are visualized in dark and light blue, respectively,  
}

\value{
  A graph showing the entropy values of the elements in pairs with top and bottom entropy scores 
}

\references{
For an application of these graphs see : 

Pele J, Abdi H, Moreau M, Thybert D and Chabbert M (2011) Multidimensional scaling reveals the main evolutionary pathways of class A G-protein-coupled receptors. \emph{PLoS ONE} \bold{6}: e19094. doi:10.1371.
}


\author{
  Julien PELE, Antoine GARNIER and Marie CHABBERT
}

\examples{
  
##Example with MSA
  #File path for output file
  wd <- tempdir()
  #wd <-getwd() 
  file <- file.path(wd,"test_seq7") 
  
  #Importing MSA file  
  msf <- system.file("msa/toy_align.msf", package = "Bios2cor")
  align <- import.msf(msf)

  #Creating OMES correlation object and selecting a correlation matrix
  correlation <- omes(align, gap_ratio = 0.2)
  corr_matrix <- correlation$Zscore 

  #Creating entropy object
  entropy <- entropy(align)

  #Creating a delta filter based on entropy
  filter <- delta_filter(entropy, Smin = 0.4, Smax = 0.7)
 
  #Creating the entropy graph
  scores_entropy.plot(entropy, corr_matrix, filepathroot = file, filter=filter)

##Example with MD
#  #File path for output file
#  wd <- tempdir()
#  #wd <-getwd() 
#  file <- file.path(wd,"test_dyn7") 

#  #Reading pdb and dcd files
#  pdb <- system.file("rotamer/toy_coordinates.pdb", package= "Bios2cor")
#  trj <- system.file("rotamer/toy_dynamics.dcd", package= "Bios2cor")

#  #Creating dynamic_structure object for selected frames 
#  wanted_frames <- seq(from = 1, to = 40, by = 1)
#  dynamic_structure <- dynamic_structure(pdb, trj, wanted_frames)
#    
#  #Creating rotamers object using conversion_file
#  conversion_file <- system.file("rotamer/dynameomics_rotamers.csv", package= "Bios2cor")
#  rotamers <- angle2rotamer(dynamic_structure, conversion_file)
# 
#  #Creating the dynamic_entropy and filter objects 
#  entropy <- dynamic_entropy(rotamers)
#  filter <- delta_filter(entropy, Smin = 0.0, Smax = 0.1)

#  #Creating correlation object
#  #dyn_cor <- dynamic_circular(dynamic_structure)
#  dyn_cor <- dynamic_omes(dynamic_structure,rotamers) 

#  #selection correlation matrix
#  corr_matrix <- dyn_cor$Zscore_noauto
# 
#  #Creating the entropy graph
#  scores_entropy.plot(entropy, corr_matrix, filepathroot = file, filter=filter)

}


