\name{bios2cor-package}
\Rdversion{2.2}
\alias{bios2cor-package}
\alias{bios2cor}
\alias{Bios2cor}
\docType{package}
\title{
Correlation/covariation Analysis in Biological Sequences and Simulations
}
\description{
The \code{bios2cor} package is dedicated to the computation and analysis of correlation/covariation between positions in multiple sequence alignments (MSA) and between side chain dihedral angles during molecular dynamics simulations (MD). Features include the computation of correlation/covariation scores using a variety of scoring functions and their analysis through a variety of tools including network representation and principal components analysis. In addition, several utility functions are based on the R graphical environment to provide friendly tools for help in data interpretation. 

For clarity purpose, version 2 of the package differentiates scoring functions working on MSA and on MD because their arguments are different. Analysis functions are common with auto detection of MSA or MD.   

  The main functionalities of \code{bios2cor} are summarized below:
  \describe{
    \item{(1) CORRELATION/COVARIATION METHODS : }{
    Methods that can be used to analyze sequence alignments and molecular simulations and to calculate a correlation/covariation matrix containing a score for each pair of positions (sequence alignment) or each pair of dihedral angles (molecular simulations).
    
    Methods working with sequence alignments (fasta or msf file is required):
      \itemize{
	\item \code{\link{omes}}: calculates the difference between the observed and expected occurrences of each possible pair of amino acids (x, y) at positions i and j of the alignment.
	
	\item \code{\link{mi}} and \code{\link{mip}}: calculate a score based on the probability of joint occurrence of events (MI)  and a corrected score by substraction of the average product (MIP), respectively. 

	\item \code{\link{elsc}}: calculates a score based on rigorous statistics of correlation/covariation in a perturbation-based algorithm. It measures how many possible subsets of size n would have the composition found in column j.

	\item \code{\link{mcbasc}}: relies on a substitution matrix giving a similarity score for each pair of amino acids.
	
      }
      
      Methods working with molecular simulations (pdb and dcd files are required) :
      \itemize{
	\item \code{\link{dynamic_circular}}: calculates a correlation/covariation score based on a circular version of the Pearson correlation coefficient, between each pair of side chain dihedral angles in a trajectory obtained from molecular dynamics simulations.
	
	\item \code{\link{dynamic_omes}}: calculates the difference between the observed and expected occurrences of each possible pair of rotamers (x, y) occuring at side chain dihedral angles i and j in a trajectory.
	
	\item \code{\link{dynamic_mi}} and \code{\link{dynamic_mip}}: calculate a score based on the probability of joint occurrence of rotameric states (MI) and a corrected score by substraction of the average product (MIP), respectively.

      }
 
      The methods working with molecular simulations require the following functions :
      \itemize{

	\item \code{\link{dynamic_structure}}: using the result of the \code{xyz2torsion} function from the \code{bio3D} package, creates a structure that contains side chain dihedral angle informations for each selected frame of the trajectory.

 	\item \code{\link{angle2rotamer}}: using the result of the \code{\link{dynamic_structure}} function, creates a structure that associates rotameric state to each side chain dihedral angle for each selected frame of the trajectory.

      }

   }
    
    \item{(2) ADDITIONNAL FUNCTIONS  : }{
    Functions that can be used to analyse the results of the correlation/covariation methods :
      
      Entropy functions :
      \itemize{
	\item \code{\link{entropy}}: calculates an entropy score for each position of the alignment. This score is based on sequence conservation and uses a formula derived from the Shannon's entropy.
	
	\item \code{\link{dynamic_entropy}}: calculates a "dynamic entropy" score for each side chain dihedral angle of a protein during molecular simulations. This score is based on the number of rotameric changes of the dihedral angle during the simulation.    
      }
      
      Filters:
      \itemize{
	\item \code{\link{delta_filter}}: given an entropy object, returns a delta filter for each position of the alignment or each side chain dihedral angle of the protein, based on entropy/dynamic entropy value.
      }
      
      PCA :
      \itemize{
	\item \code{\link{centered_pca}}: returns a principal component analysis of the double-centered correlation/covariation matrix passed as a parameter. A delta filter can be precised.
      }
    }
    
    \item{(3) OUTPUT FILES : }{
      Functions that can be used to produce output files.
       
      Some data structures can be stored in txt/csv files :
      \itemize{
	\item \code{\link{write.scores}}: Using the result of a correlation/covariation method, creates a file containing the score of each pair of positions (sequence alignment analysis) or of side chain dihedral angles (molecular simulations) and optionaly their entropy/dynamic_entropy score.

	\item \code{\link{top_pairs_analysis}}: Using the result of a correlation/covariation method and an integer N, creates two files containing (1) the top N pairs with their scores and (2) the individual elements, their contact counts and their entropy score for the top N pairs. Subsequently, these files can be  used for network visualization with the Cytoscape program accessible at \emph{https://cytoscape.org}. 
 
	\item \code{\link{write.pca}}: Using the result of the \code{\link{centered_pca}} function, creates a file that contains the coordinates of each element in the principal component space.
 
	\item \code{\link{write.pca.pdb}}: Using the result of the \code{\link{centered_pca}} function, creates a pdb file with the PCA coordinates on three principal components along with a pml file for nice visualization with the Pymol molecular visualization program accessible at \emph{https://pymol.org}.
      }
      
      Some data can be visualized in png/pdf files:
      \itemize{
	\item \code{\link{scores.boxplot}}: Using the result of one or several correlation/covariation methods, creates a boxplot to visualize the distribution of the Z-scores.

	\item \code{\link{network.plot}}: Using the result of the \code{\link{top_pairs_analysis}} function, creates the graph of a network representation of the data.

	\item \code{\link{scores_entropy.plot}}: Using the result of a correlation/covariation method and an entropy structure, creates a graph comparing correlation/covariation scores with entropy values. Each pair of elements (i,j) is placed in the graph with (entropy[i] ; entropy[j]) as coordinates. The color code of each point is based on its correlation/covariation score (red/pink color for top values, blue/skyblue for bottom values).

	\item \code{\link{pca_screeplot}}: Using the result of the \code{\link{centered_pca}} function, creates the graph of the eigen values (positive values only).

	\item \code{\link{pca_2d}}: Using the result of the \code{\link{centered_pca}} function, creates a graph with the projection of the elements on two selected components.

	\item \code{\link{angles.plot}}: Using pdb and dcd files and the result of a correlation/covariation method, creates graphs to monitor the time evolution of each dihedral angle in the top N pairs
      }
      
    }
  }
}
\details{
\tabular{ll}{
  Package: \tab BioCor\cr
  Type: \tab Package\cr
  Version: \tab 2.1\cr
  Date: \tab 2020-01-30\cr
  License: \tab GPL \cr
}
}
\author{
  Bruck TADDESE [aut], Antoine Garnier [aut], Madeline DENIAUD [aut], Lea BELLANGER [ctb], Julien PELE[ctb], Jean-Michel BECU [ctb], Marie CHABBERT [cre].
  Maintainer: Marie CHABBERT <marie.chabbert@univ-angers.fr>
}


\examples{
  #File path for output files
  wd <- tempdir()
  #wd <- getwd()
  file <- file.path(wd,"test_seq1") 

  #Importing multiple sequence alignment
  msf <- system.file("msa/toy_align.msf", package = "Bios2cor")
  align <- import.msf(msf)

  #Creating correlation object with OMES method
  omes <- omes(align, gap_ratio = 0.2)

  #Creating entropy object
  entropy <- entropy(align, gap_ratio = 0.2)

  #Creating delta filter based on entropy
  filter <- delta_filter(entropy, Smin = 0.3, Smax = 0.8)
  
  #Selecting a correlation matrix
  omes <-omes$score

  # Creating PCA object for selected correlation matrix and storing eigen values in csv file
  pca <- centered_pca(omes, filepathroot= file, pc= NULL, dec_val= 5, filter= filter)
  
  
}
