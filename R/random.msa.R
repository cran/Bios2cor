# Package: Bios2cor 
# This file is part of the Bios2cor and Bios2mds R package.
# Bios2cor and Bios2mds are free softwares: you can redistribute them and/or modify
# them under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# See the GNU General Public License at:
# http://www.gnu.org/licenses/
#

random.msa <- function (nb.seq = 100, id = "SEQ", nb.pos = 100, gap = FALSE, aa.strict = FALSE,align = NULL, align.replace = TRUE) {

  #one letter codes for amino acids 
  aa <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K",
    "M", "F", "P", "S", "T", "W", "Y", "V", "B", "Z", "J", "X")
  replace <- TRUE
  #remove ambiguous amino acids
  if (aa.strict)
    aa <- aa[1:20]

  if (gap)
    aa <- c(aa, "-")

  if(!is.null(align)){
    if (!inherits(align, "align"))
        stop("mmds is not a 'align' object")
    aa <-as.vector(unlist(align))
    replace<-align.replace
  }

  msa <- lapply(seq_len(nb.seq), function (i) {sample(aa, nb.pos, replace = replace)})

  msa.names <- paste(id, seq_len(nb.seq), sep = "")

  names(msa) <- msa.names

  return(msa)
}
