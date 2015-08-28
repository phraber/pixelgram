## pixgramr
####An R package that implements paired tree (phylogram) and pixel plots to view polymorphisms among aligned sequences.

####Synopsis
This package will plot a nucleotide and/or protein sequence alignment together with a phylogenetic tree.  Alignments are shown in pixel (raster) format, whereby sequences correspond to rows and follow the tree leaf order.  Options are provided to specify a reference sequence for numbering positions and master sequence (or consensus) to transform the pixel view, whereby only sites that mismatch the master sequence are shown.  This enables viewing sequence variation in light of the phylogeny and vice versa.  Many other options exist to decorate the tree, which is rendered using the ape package.

####Quick Start
Any R user can install this package with just a few lines lines of code:
1. Start R, by whatever method you choose.
1. Type "install.packages('devtools')"
1. If prompted to compile a required package from source, you can safely decline without ill effects.
1. Type "devtools::install_github('phraber/pixgram')".
1. Type "vignette('pixgramr')"

If you cannot view the vignette, please see this PDF file for examples - https://github.com/phraber/pixgram/blob/master/examples.pdf

#### To list all functions provided in this package:
1. Start R, by whatever method you choose.
1. Type "library(pixgramr)".
1. Type "ls('package:pixgram')".
