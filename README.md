## pixgramr
####An R package that implements paired tree (phylogram) and pixel plots to view polymorphisms among aligned sequences.

####Synopsis
This package will plot a nucleotide and/or protein sequence alignment together with a phylogenetic tree.  Alignments are shown in pixel (raster) format, whereby sequences correspond to rows and follow the tree leaf order.  Options are provided to specify a reference sequence for numbering positions and master sequence (or consensus) to transform the pixel view, whereby only sites that mismatch the master sequence are shown.  This enables viewing sequence variation in light of the phylogeny and vice versa.  Many other options exist to decorate the tree, which is rendered using the ape package.

####Quick Start
Any R user can install this package with just a few lines lines of code:

  install.packages("devtools")

  devtools::install_github("phraber/pixgram")

  library(pixgramr)
