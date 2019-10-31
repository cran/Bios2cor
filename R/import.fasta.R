import.fasta <- function (file, aa.to.upper = TRUE, gap.to.dash = TRUE,log.file=NULL) {

  if(missing(file)) {
    stop("file is missing")
  }

  # Read as a vector of lines
  lines <- readLines(file)

  # Localize sequence identifiers and check fasta format
  loc <- grep(">", lines)
  if (length(loc) == 0){
	if(!is.null(log.file))
		write("file is not in fasta format",log.file)  
        stop("file is not in fasta format")
  }
  # Get sequence identifiers
  id <- sub("^>(\\S+).*$","\\1", lines[loc])
  nb.seq <- length(id)

  # Localize sequence pieces for each identifier
  start <- loc + 1
  end <- loc - 1  
  end <- c(end[-1], length(lines))

  seq <- sapply(seq_len(nb.seq), function(i) {paste(lines[start[i]:end[i]], collapse = "")})
  seq <- gsub("\\s", "", seq)

  # Turn aa into upper case
  if (aa.to.upper)
    seq <- toupper(seq)

  # Give a list of split sequences
  seq <- strsplit(seq, split = "")
  names(seq) <- id

  # Turn gap into dash character
  if (gap.to.dash)
    seq <- lapply(seq, function (i) {i[is.gap(i)] <- "-"; return(i)})
  class (seq) <- c("align")
  return(seq)
}
