\name{bios2cor-package}
\Rdversion{1.1}
\alias{bios2cor-package}
\alias{bios2cor}
\alias{Bios2cor}
\docType{package}
\title{
Correlation Analysis in Biological Sequences and Simulations
}
\description{
Bios2cor is dedicated to the computation and analysis of correlation/co-variation between positions in multiple sequence alignments and between side chain dihedral angles during molecular dynamics simulations.  Features include the ability to compute correlation/co-variation using a variety of scoring functions and to analyze the correlation/co-variation matrix through a variety of tools including network representation and principal components analysis.  In addition, several utility functions are based on the R graphical environment to provide friendly tools for help in data interpretation. 

  The main functionalities of \code{Bios2cor} are summarized below:
  \describe{
    \item{(1) CORRELATION/COVARIATION METHODS : }{
    Methods that can be used to analyze sequence alignments nad molecular simulations and to calculate a correlation matrix containing a score for each pair of positions (sequence alignment) or each pair of dihedral angles (molecular simulations).
    
    Methods working with sequence alignments (fasta or msf file is required):
      \itemize{
	\item \code{\link{omes}}: calculates the difference between the observed and expected occurrences of each possible pair of amino acids (x, y) at positions i and j of the alignment.
	
	\item \code{\link{mip}}: calculates a score (MI) based on the probability of joint occurrence of events and correct it with the average product correction which is subtracted from the MI score.

	\item \code{\link{elsc}}: calculates a score based on rigorous statistics of covariation in a perturbation-based algorithm. It measures how many possible subsets of size n would have the composition found in column j.

	\item \code{\link{mcbasc}}: relies on a substitution matrix giving a similarity score for each pair of amino acids.
	
      }
      
      Methods working with molecular simulations (pdb and dcd files are required) :
      \itemize{
	\item \code{\link{rotamer_circular}}: calculates a correlation score based on a circular version of the Pearson correlation coefficient, between each pair of side chain dihedral angles in a trajectory obtained from molecular dynamics simulations
	
	\item \code{\link{rotamer_omes}}: calculates the difference between the observed and expected occurrences of each possible pair of rotamers (x, y) occuring at side chain dihedral angles i and j in a trajectory
	
	\item \code{\link{rotamer_mip}}: calculates a score (MI) based on the probability of joint occurrence of rotameric states and correct it with the average product which is subtracted from this MI score.

      }
 
      The methods working with molecular simulations require the following functions :
      \itemize{

	\item \code{\link{dynamic_struct}}: using the result of the \code{xyz2torsion} function, creates a unique structure that contains side chain dihedral angle informations for each selected frame of the trajectory

 	\item \code{\link{angle_conversion}}: using the result of the \code{\link{dynamic_struct}} function, creates a structure that associates rotameric state to each side chain dihedral angle for each selected frame of the trajectory.

      }

   }
    
    \item{(2) ADDITIONNAL FUNCTIONS  : }{
    Functions that can be used to analyse the results of the correlation methods :
      
      Entropy functions :
      \itemize{
	\item \code{\link{entropy}}: calculates an entropy score for each position of the alignment. This score is based on sequence conservation and uses a formula derived from the Shannon's entropy.
	
	\item \code{\link{rotamer_entropy}}: calculates a "dynamic entropy" score for each side chain dihedral angle of a protein during molecular simulations. This score is based on the number of rotameric changes of the dihedral angle during the simulation. 
      }
      
      Filters:
      \itemize{
	\item \code{\link{gauss_weighting}}: given an entropy object, returns a \emph{"gaussian weighting"} for each position of the alignment or each side chain dihedral angle of the protein.
	\item \code{\link{sigmoid_weighting}}: given an entropy object, returns a \emph{"sigmoidal weighting"} for each position of the alignment or each side chain dihedral angle of the protein.
	\item \code{\link{delta_weighting}}: given an entropy object, returns a \emph{"delta weighting"} for each position of the alignment or each side chain dihedral angle of the protein.
      }
      
      PCA :
      \itemize{
	\item \code{\link{centered_pca}}: returns a principal component analysis of the correlation matrix passed as a parameter. A weighting filter can be precised, using the weighting functions described above.
      }
    }
    
    \item{(3) OUTPUT FILES : }{
      Functions that can be used to produce output files (txt/csv and graphs)
      
      Some data structures can be stored in txt/csv files :
      \itemize{
	\item \code{\link{create_corrfile}}: Using the result of a correlation/covariation method, creates a file containing the score of each pair of positions (sequence alignment analysis) or of side chain dihedral angles (molecular simulations)
	\item \code{\link{corr_contact}}: Using the result of a correlation/covariation method and an integer X, creates a file containing top positions/dihedral angles and their contact counts (number of times the position/dihedral angle appears in the top X pairs with highest scores 
	\item \code{\link{create_pcafile}}: Using the result of the \code{\link{centered_pca}} function, creates a file that contains the coordinates of each position or of each dihedral angle in the principal components
	\item \code{\link{write_pdb}}: Using the result of the \code{\link{centered_pca}} function, creates a pdb file with the PCA coordinates on three principal components along with a pml file for nice visualization with Pymol
      }
      
      Some data can be analyzed thanks to special graphs :
      \itemize{
	\item \code{\link{create_boxplot}}: Using the result of a correlation/covariation method, creates a boxplot to visualize the distribution of the Z-scores
	\item \code{\link{create_network}}: Using the result of the \code{\link{corr_contact}} function, creates the graph of a network representation of the data with links between correlated/covarying elements
	\item \code{\link{entropy_graph}}: Using the result of a correlation/covariation method and an entropy structure, creates a graph comparing correlation scores with entropy values. Each pair of elements (i,j) is placed in the graph with (entropy[i] ; entropy[j]) as coordinates. The color code of each point is based on its correlation score (red/pink color for top values, blue/skyblue for bottom values).
	\item \code{\link{create_screeplot}}: Using the result of the \code{\link{centered_pca}} function, creates the graph of the eigen values (positive values only)
	\item \code{\link{pca_2d}}: Using the result of the \code{\link{centered_pca}} function, creates a graph with the projection of the elements on two selected components
	\item \code{\link{ang_evo_graph}}: Using pdb and dcd files and the result of a correlation/covariation method, creates graphs to monitor the time evolution of each dihedral angle in the top X pairs
      }
      
    }
  }
}
\details{
\tabular{ll}{
  Package: \tab BioCor\cr
  Type: \tab Package\cr
  Version: \tab 1.5\cr
  Date: \tab 2017-03-06\cr
  License: \tab GPL \cr
}
}
\author{
  Bruck TADDESE [aut], Antoine Garnier [aut], Madeline DENIAUD [aut], Lea BELLANGER [ctb], Julien PELE[ctb], Jean-Michel BECU [ctb], Marie CHABBERT [cre].
  Maintainer: Marie CHABBERT <marie.chabbert@univ-angers.fr>
}

\keyword{correlation}
\keyword{methods}
\keyword{alignment}

\examples{
  msf <- system.file("msa/toy_align.msf", package = "Bios2cor")
  align <- import.msf(msf)

  #Creating OMES object
  omes <- omes(align,fileHelix= NULL , diag= 0, fileCSV= NULL, gap_val= 0.8, z_score= TRUE)

  #Creating ENTROPY object
  entropy <- entropy(align)

  #Creating weighting filter
  filter <- gauss_weighting(entropy, L= 0.1)
  
  # Creating PCA structures for OMES method and storing in txt file
  omes <-omes$normalized
  pca <- centered_pca(omes, m= filter, pc= NULL, dec_val= 5,eigenvalues_csv= NULL)
  
  
}
