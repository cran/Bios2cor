\name{entropy_graph}

\Rdversion{1.1}

\alias{entropy_graph}

\title{
 Entropy graph
}

\description{
  The aim of this graph is the visualization of the entropy values of the elements in pairs with top and bottom entropy scores 
}

\usage{
  entropy_graph(entropy, corr, filter= NULL, elite, high, csv, name)
}

\arguments{
 \item{entropy}{
  An object created by the \bold{\emph{entropy}} (sequence alignment) or the \bold{\emph{rotamer_entropy}} (trajectory) function. 
 }
 \item{corr}{
  An object created by one of the correlation functions: 
 }
 \item{filter}{
  A vector created by one of the weighting functions (gaussian, delta or sigmoid filter)
 }
 \item{elite}{
  An integer to determine the number of pairs with the highest and lowest scores (e.g. 25: pairs ranked 1 to 25 in decreasing or increasing order) to be colored with the "elite" color codes. 
 }
 \item{high}{
  An integer to determine the number of pairs with the next highest and lowest scores (e.g. 275: pairs ranked 26 to 300 in decreasing or increasing order) to be colored with the "high" color codes.
 }
 \item{csv}{
  Logical value that indicates if information must be stored in a csv file
 }
 \item{name}{
  The full path name of the graph that will be created
 }
}

\details{
  Using the result of a correlation/covariation method and an entropy structure, creates a graph comparing correlation scores with entropy values. Each pair of elements (i,j) is placed in the graph with (entropy[i] ; entropy[j]) as coordinates. The color code of each point is based on its correlation score (blue/red color for elite values, light blue/pink for bottom values).
}

\value{
  A graph showing the entropy values of the elements in pairs with top and bottom entropy scores 
}

\author{
  Julien PELE and Antoine GARNIER
}

\examples{
  msf <- system.file("msa/toy_align.msf", package = "Bios2cor")
  align <- import.msf(msf)

  #Creating OMES object
  omes <- omes(align,fileHelix= NULL , diag= 0, fileCSV= NULL, gap_val= 0.8, z_score= TRUE)
  omes <- omes$normalized
  
  #Creating ENTROPY object
  entropy <- entropy(align)
    
  #creating entropy graph
  entropy_graph(entropy, omes, filter = NULL, elite = 25, high = 275, csv = TRUE, 
          name = "entropy_graph.pdf")
  
}

\keyword{entropy}
\keyword{graph}
