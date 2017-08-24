\name{create_boxplot}

\Rdversion{1.1}

\alias{create_boxplot}

\title{
  Creates boxplot file 
}

\description{
  Given correlation objects, creates a boxplot file to compare themselves
}

\usage{
  create_boxplot(correlation_objects, objects_names, elite, high, filepath)
}

\arguments{
 \item{correlation_objects}{
  Correlation objects to be compared
 }
 \item{objects_names}{
  Correlation objects names
 }
 \item{elite}{
  Integer that will define the number of top and bottom values of correlation objects to be taken into account (e.g. top and bottom 25)
 }
 \item{high}{
  Integer that will define the number of next top and bottom values of correlation objects to be taken into account (e.g. top and bottom 275)
 }
 \item{filepath}{
  The full path name of the output file 
 }
}

\details{
  Correlation structures contain correlation/covariation scores for each pair of elements [i,j].
  The boxplot will compare these structures using color codes : the highest values are dark blue, the next highest values are light blue, the lowest values are red and the next lowest values are pink.
}

\value{
  A boxplot comparing different correlation objects
}

\author{
  Julien PELE and Antoine GARNIER
}

\examples{

msf <- system.file("msa/toy_align.msf", package = "Bios2cor")
  align <- import.msf(msf)

  #Creating OMES object
  omes <- omes(align,fileHelix= NULL , diag= 0, fileCSV= NULL, gap_val= 0.8, z_score= TRUE)
  omes <-omes$normalized

  #Creating list for boxplots
    corr_objects <- list(omes)
    corr_objects_names <- c("omes")
    create_boxplot(corr_objects, corr_objects_names, 25, 275,"boxplot.png")
}

\keyword{boxplot}
\keyword{correlation}
