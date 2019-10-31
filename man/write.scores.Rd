\name{write.scores}

\Rdversion{1.1}

\alias{write.scores}

\title{
  Creates a correlation/covariation ouput file 
}

\description{
  Given a correlation object, writes the score, the Zscore and, optionally the entropy, for each pair of elements in a csv file.
}

\usage{
  write.scores(correlation, entropy= NULL, filepathroot=NULL)
}

\arguments{
 \item{correlation}{
  An object created by a correlation/covariation function (\code{\link{omes}}, \code{\link{mi}}, \code{\link{mip}}, \code{\link{elsc}}, \code{\link{mcbasc}}, \code{\link{dynamic_circular}}, \code{\link{dynamic_omes}}, \code{\link{dynamic_mi}}, \code{\link{dynamic_mip}}
 }
 \item{entropy}{
  An object created by the \code{\link{entropy}} function 
 }
 \item{filepathroot}{
  The root of the full path name for the output file. DEfault is NULL (a "CORR_SCORES.csv" file is created). If not NULL, the "_CORR_SCORES.csv" extension is added to the root name.
 }
}

\details{
  Elements represent positions in multiple sequence alignments and side chain dihedral angles in molecular dynamic simulations (MD). 

  In sequence analysis, the correlation object contains two matrices with the correlation/covariation scores and Z-scores for each pair of elements [i,j].   
  If entropy = NULL, each line of the output file will contain element  i, element  j, score[i,j], and Z-score[i,j]. 
  If entropy is not NULL, each line of the output file will contain element  i, element  j, score[i,j], Zscore[i,j], entropy[i], and entropy[j].  

  In MD analysis, the correlation object contains four matrices with (1) the correlation/covariation scores for each pair of rotamers (score), (2) the Z-scores for each pair of rotamers (Zscore), (3) the correlation/covariation scores for each pair of rotamers with zero values for autocorrelation (correlation within the same side chain) (score_noauto) and (4) the Z-scores calculated without autocorrelation pairs and zero values for autocorrelation pairs (Zscore_noauto). 
  If entropy = NULL, each line of the output file will contain element  i, element  j, score[i,j], Zscore[i,j], score_noauto[i,j], and Zscore_noauto[i,j].
  If entropy is not NULL, each line of the output file will contain element  i, element  j, score[i,j], Zscore[i,j], score_noauto[i,j], Zscore_noauto[i,j], entropy[i], and entropy[j].

}


\value{
  A csv file containing the correlation/covariation scores and optionally the entropies.
}

\author{
  Antoine GARNIER and Marie CHABBERT
}

\examples{

#Example for MSA 
  #File path for output files
  wd <- tempdir()
  #wd <-getwd() 
  file <- file.path(wd,"test_seq") 

  #Importing MSA file
  msf <- system.file("msa/toy_align.msf", package = "Bios2cor")
  align <- import.msf(msf)

  #Creating correlation and entropy objects
  correlation <- omes(align, gap_ratio= 0.2)
  entropy <- entropy(align)
  
  #Writing results to csv file
  write.scores(correlation, entropy, file)

###Example for MD
  #File path for output files
  wd <- tempdir()
  #wd <-getwd() 
  file <- file.path(wd,"test_dyn") 

  #Reading the pdb and dcd files and the angles to rotamers conversion file
  pdb <- system.file("rotamer/toy_coordinates.pdb", package= "Bios2cor")
  trj <- system.file("rotamer/toy_dynamics.dcd", package= "Bios2cor")
  conversion_file <- system.file("rotamer/dynameomics_rotamers.csv", package= "Bios2cor")

  #Creating dynamic_structure and rotamers objects
  wanted_frames <- seq(from = 5, to = 40, by = 10)
  dynamic_structure <- dynamic_structure(pdb, trj, wanted_frames)
  rotamers <- angle2rotamer(dynamic_structure, conversion_file)
  
  #Creating correlation and entropy objects
  wanted_residues <- c("F","H","N","Y","W")
  #dyn_corr <- dynamic_omes(dynamic_structure, rotamers, wanted_residues)
  dyn_corr <- dynamic_circular(dynamic_structure, wanted_residues)
  dyn_entropy <- dynamic_entropy(rotamers)
   
  #Writing results to csv file
  write.scores(dyn_corr, dyn_entropy, file)

}


