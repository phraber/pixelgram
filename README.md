# pixgram
###An R package that implements paired tree (phylogram) and pixel plots to view polymorphisms among aligned sequences.

####Synopsis
This package will plot a nucleotide and/or protein sequence alignment together with a phylogenetic tree.  Alignments are shown in pixel (raster) format, whereby sequences correspond to rows and follow the tree leaf order.  Options are provided to specify a reference sequence or consensus to transform the pixel view, whereby only sites that mismatch the reference sequence are shown.  This enables viewing sequence variation in light of the phylogeny.  Many other options exist to decorate the tree, which is rendered using the ape package.
