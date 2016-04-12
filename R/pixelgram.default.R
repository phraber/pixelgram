#' @export
pixelgram.default <- function(...) UseMethod("pixelgram")

#' Create (construct) a pixelgram object.
#'
#' \code{pixelgram()} returns a pixelgram object for plotting.
#'
#' At minimum, a nucleotide sequence alignment must be provided. If no tree is given, \code{ape::bionj} will make a neighbor-joining tree with simplified parameters.  The alignment can be translated (in +1 reading frame) with the standard genetic code by setting \code{is_orf=T}.  For other translation options see \code{seqinr::translate}.
#'
#' The \code{refseq_name} argument is used to specify standardized sequence numbering.  By default it is HXB2, which will be excised from the alignment after creating a lookup-table to number positions in the alignmnent.  If the exact string match to this argument is not found among sequence names in the alignment, it is searched as a regular expression, but the resulting hit is not excised from the alignment.  (This is illustrated in the example vignette, in which HXB2 with a long name is not excised from the subtype reference alignment.)  If \code{excise_refseq=F} the reference sequence is left in the alignment and could also become the master sequence.
#'
#' The \code{master_name} argument is used when applying transformations to the pixel plot.  These transformations follow the behavior of highlighter and pixel plot utilities at the LANL HIV database, i.e. only differences from the master sequence are shown in each column.  The default is to use the first non-reference sequence in the alignment as the master.  That is, the master sequence is set after having excised the reference sequence.  In addition to names in the sequence alignment, alternative options for the master sequence include "__none__", "__first__" (the default), "__last__", and "__consensus__".
#'
#' The option \code{xform_type} is an integer from 0 to 3.  The default (0) means no transformation is applied, yielding a pixel plot that matches the "default" color scheme in the online pixel utility.  That is, positions in all rows are colored.  Other transformations will change matches to each master site as background "transparent" colors or mutational changes using an extensible color scheme.  For details about how the color scheme is defined, run \code{help(pixmap_colors)}.  Option 3 works for amino-acid alignments to show changes in charge and potentially N-linked glycosylation sites, relative to the master sequence, but is slow to compute.  Setting the option \code{xform_master=T} applies the transformation to the master sequence (following the Highlighter plot scheme, the default is F).
#'
#' @param tre_file File containing Newick-formatted tree.
#' @param tre ape phylo object in memory
#' @param nts_file File containing nucleotide sequence alignment.
#' @param nts seqinr nucleotide alignment in memory.
#' @param aas_file File containing protein sequence alignment.
#' @param aas seqinr protein alignment in memory.
#' @param alignment_format Format of alignment file/s; must be one of these: \code{"fasta"}, \code{"clustal"}, \code{"phylip"}, \code{"msf"}, or \code{"mase"}.
#' @param is_orf flag to indicate the nt alignment is +1 codon aligned
#' @param pngs2o Switch to mark asparagines (N) in PNG motifs as O.
#' @param refseq_lut Look-up-table for reference numbering of alignment positions
#' @param refseq_name Reference sequence name in the alignment.  (Do not confuse "reference" sequence, used for numbering, with "master" sequence, used for coloring.)
#' @param excise_refseq Should the reference sequence be removed from the alignment before rendering?
#' @param master_name Name of sequence to use when showing only polymorphisms, i.e. sites that do not match. Must be a name in the alignment or one of \code{"__first__"}, \code{"__last__"}, \code{"__consensus__"}, or \code{"__none__"}.  (Try not to confuse "master" and "reference" sequences.  They could be the same sequence, but are used differently.)
#' @param xform_master Should the master sequence be transformed when rendered in the pixel plot?
#' @param xform_type Controls the way polymorphisms are shown in pixel plot; see Details.
#' @param show_tree Show the tree?  If not, only show the pixel plot.
#' @param raster_width Proportionate width of the pixel plot, e.g. 1 indicates equal width.
#' @param raster_margin Offset in raster portion of the layout, e.g. for heatmap.
#' @param x_lim Range of x coordinates to use for plotting, set automatically.
#' @param y_lim Range of y coordinates to use for plotting, set automatically.
#' @param invert_y Should the y-axis be inverted?
#' @param main Plot title, e.g. coded subject id.  If given, it is plotted at the top of the pixelgram output.
#' @param sub Plot subtitle, i.e. region sequenced.  If specified, it appears with "site" below the pixel plot.
#'
#' @return pixelgram object
#'
#' @seealso
#'
#' External methods used to translate nts to aas and build a tree from nts:
#' \code{\link{seqinr}{translate}}, \code{\link{ape}{bionj}}, and \code{\link{ape}{dist.dna}}
#'
#' Online examples that are precursors to this package:
#' \url{http://hiv.lanl.gov/content/sequence/pixel/pixel.html} and
#' \url{http://hiv.lanl.gov/content/sequence/HIGHLIGHT/highlighter_top.html}.
#'
#' @examples
#' \dontrun{
#'   plot( pg <- pixelgram::pixelgram(nts=pixelgram::hiv.ref$nts) )
#' }
#'
#' @export

pixelgram <- function(tre_file=NULL,
                            tre=NULL,
                            nts_file=NULL,
                            nts=NULL,
                            aas_file=NULL,
                            aas=NULL,
			    alignment_format="fasta",
			    is_orf=F,
			    pngs2o=F,

                            refseq_lut=NULL,
			    refseq_name="HXB2",
			    excise_refseq=T,

                            master_name="__first__", # was __none__
                            xform_master=F,

                            xform_type=0,
                            show_tree=T,
#                            palette=NULL, # used for tree
                            raster_width=1,
                            raster_margin=0.025,

			    x_lim=NULL,
			    y_lim=NULL,

                            invert_y=F,
main=NULL, sub=NULL) {

    message("*** Creating pixelgram ***")

    P <- list(tre_file=tre_file,
              nts_file=nts_file,
              aas_file=aas_file,
	      alignment_format=alignment_format,
              is_orf=is_orf,
	      pngs2o=pngs2o,
              tre=tre,
              nts=nts,
              aas=aas,

              refseq_lut=refseq_lut,
	      refseq_name=refseq_name,
	      excise_refseq=excise_refseq,

              master_name=master_name,
#              master_index=master_index,
              nt_master=NULL, #vector that contains the master nt sequence
              aa_master=NULL, #vector that contains the master aa sequence
              xform_master=xform_master,

#	      pixmap_colors=pixmap_colors, # this is used in the pixel plot
              nt_rast=NULL,#nts,
              aa_rast=NULL,#aas,

              x_lim=x_lim,
              y_lim=y_lim,
              invert_y=invert_y,
###              cex=NULL,

 #             nt_sites=NULL,
 #             aa_sites=NULL,
              xform_type=xform_type,
              show_tree=show_tree,
              raster_width=raster_width,
              raster_margin=raster_margin,
              main=main, sub=sub)

              # prefix?
              #	main=main,
              #	sub=sub,
#              strip_columns=NULL, #strip_columns,
#              palette=NULL)#, #palette, # this is used for tree decoration
    #        verbose=verbose)

    class(P) <- "pixelgram"

    if (!is.null(tre_file))
        P <- set.tre.file(P, tre_file)

    if (!is.null(nts_file))
        P <- set.aln.file(P, nts_file,
            ifelse(P$is_orf, "codon", "nt"), P$alignment_format)

    if (!is.null(P$nts))
	P <- set.refseq(P, is.aa=F) # or it may never be invoked

    if (P$is_orf & !is.null(P$nts))
	P <- translate.codons(P)

    if (!is.null(aas_file))
        P <- set.aln.file(P, aas_file, "aa", P$alignment_format)

    if (!is.null(P$aas))
	P <- set.refseq(P, is.aa=T) # otherwise this may never get called

    return ( P )  ### BEWARE that this has not yet been validated
}

#    if (strip_columns) {
#	P <- pixelgram.colstrip(P)
#    } else {
#        # could help later with refseq lookups?
#        if (!is.null(P$aas))
#            P$aa_sites = c(1:ncol(aas))
#        if (!is.null(P$nts))
#            P$nt_sites = c(1:ncol(nts))
#    }

#' Sets the alignment file, whether nt or aa.  Also sets the refseq.
#'
#' @param P pixelgram object
#' @param f alignment file name
#' @param file_type Specifies alignment file type.  Must be "aa", "nt", or "codon".
#' @param alignment_format Alignment file format, one "fasta", "clustal", "phylip", "msf", or "mase"
#' @param pngs2o Switch to mark asparagines (N) in PNG motifs as O.
#' @return updated pixelgram object
#' @export
set.aln.file <- function(P, f, file_type='aa', alignment_format='fasta', 
    pngs2o=NULL) {

    if (class(P) != "pixelgram")
	stop("set.aln.file ERROR: Please specify a valid pixelgram object")

    if (! file_type %in% c("nt", "aa", "codon"))
	stop("set.aln.file ERROR: Please specify a valid file type")

    if (file_type != "aa") {

	P$nts_file = f

        P = set.nt_aln.from.file(P, alignment_format)

        if (file_type == "codon" | P$is_orf)
	    P <- translate.codons(P)

	P = set.refseq(P, is.aa=F)
#	P$nt_rast = P$nts

#	if (file_type == "codon") # set after refseq excision not in the conditional above
#	    P$aa_rast = P$aas

    } else {
	P$aas_file = f
        P = set.aa_aln.from.file(P, alignment_format)
	P = set.refseq(P, is.aa=T)
#	P$aa_rast = P$aas
    }

    if (!is.null(pngs2o))
	P$pngs2o = pngs2o

    if (P$pngs2o)
	P$aas <- pngs2o(P$aas)
    else
	P$aas <- unset.pngs2o(P$aas)

    return ( P )
}

#' @keywords internal
set.aa_aln.from.file <- function(P, alignment_format='fasta') {

    if (!alignment_format %in% c('fasta', 'clustal', 'phylip', 'msf', 'mase'))
        stop("set.aa_aln.from.file WARNING: specify a valid aa alignment file format")

    if (!is.null(P$aas))
        warning("set.aa_aln.from.file WARNING: updating existing aa alignment")

    if (!is.null(P$aas_file)) {

	if (file.exists(P$aas_file))
	    try(aas_seqinr <- seqinr::read.alignment(P$aas_file, alignment_format))
        else 
	    stop(paste("pixelgram ERROR: Cannot find file", P$aas_file))

        if (!is.null(aas_seqinr))
            P$aas <- tolower(seqinr::as.matrix.alignment(aas_seqinr))

	if (!is.null(P$aas))
	    message(paste("Read aa alignment file", P$aas_file))
    }
### note that if the file is null we simply pass through and keep running
    return ( P )
}

#' @keywords internal
set.nt_aln.from.file <- function(P, alignment_format='fasta') {

    if (!alignment_format %in% c('fasta', 'clustal', 'phylip', 'msf', 'mase'))
        stop("set.nt_aln.from.file WARNING: specify a valid aa alignment file format")

    if (!is.null(P$nts))
        warning("set.nt_aln.from.file WARNING: updating existing nt alignment")

    if (!is.null(P$nts_file)) {

	if (file.exists(P$nts_file))
	   try(nts_seqinr <- seqinr::read.alignment(P$nts_file, alignment_format))

        if (!is.null(nts_seqinr))
            P$nts <- tolower(seqinr::as.matrix.alignment(nts_seqinr))

	if (!is.null(P$nts))
            message(paste("Read nt alignment file", P$nts_file))
    }
### note that if the file is null we simply pass through and keep running
    return ( P )
}

#' Assign a Newick-formatted tree file to a pixelgram object.
#'
#' @param x Pixelgram object
#' @param f Name of file containing a Newick tree.
#' @return Updated pixelgram object
#' @export
set.tre.file <- function(x, f) {

    if (class(x) != "pixelgram")
	stop("set.tre.file ERROR: Please specify a valid pixelgram object")

    if (!is.null(f))
        x$tre_file = f

    x = set.tree.from.file(x)

    return ( x )
}

#' @keywords internal
translate.codons <- function(P) {

    if (class(P) != "pixelgram")
	stop("translate.codons ERROR: Please specify a valid pixelgram object")

    if (!is.null(P$nts) & is.null(P$aas)) {

	if (ncol(P$nts) %% 3 != 0)
	    warning("Attempting to translate although the number of nucleotides is not a multiple of three.")

        message("    translating")

        P$aas <- matrix(NA, ncol=floor(ncol(P$nts)/3), nrow=nrow(P$nts))

        for (i in 1:nrow(P$nts))
            P$aas[i, ] <- tolower(seqinr::translate(P$nts[i, ], ambiguous=F))
            # TODO: enable option to specify alternative/non-standard translation tables via numcode=

#            P$aas[i, ] <- tolower(seqinr::getTrans(P$nts[i, ], ambiguous=T))
	rownames(P$aas) <- rownames(P$nts)
#        colnames(P$aas) <- c(1:ncol(P$aas))

	P$is_orf = T
    }

    return ( P )
}

#' @keywords internal
set.tree.from.file <- function(P) {

    if (class(P) != "pixelgram")
	stop("set.tree.from.file ERROR: Please specify a valid pixelgram object")

    if (!is.null(P$tre))
        warning("set.tree.from.file WARNING: updating existing tree")

    if (!is.null(P$tre_file))
        if (file.exists(P$tre_file))
	    try(P$tre <- ape::read.tree(P$tre_file))
## again note that we pass through intact if the file is null
    return ( P )
}

#' @keywords internal
pixelgram.colstrip <- function(P) {

    if (class(P) != "pixelgram")
	stop("pixelgram.colstrip ERROR: Please specify a valid pixelgram object")

    message("*** Stripping ***")

    R <- P

    class(R) <- "pixelgram"

    if (!is.null(P$aas)) {

	aa_colsum <- sapply(1:ncol(P$aas), function(i)
	    length(table(P$aas[, i])))

	aa_sites <- which(aa_colsum > 1) # 2 min states may not work for all xformations

	R$aas = P$aas[ , aa_sites]
#	R$aa_rast = P$aa_rast[ , aa_sites]
	R$aa_cols = aa_sites
    }

    if (!is.null(P$nts)) {

	nt_colsum <- sapply(1:ncol(P$nts), function(i)
	    length(table(P$nts[, i])))

	if (!is.null(P$aas)) {

	    # need to maintain intact codons
	    nt_keep_sites <- rep(F, length(nt_colsum))

	    for (i in 1:length(aa_colsum)) {

		nt_position = ((i - 1) * 3) + 1

		if (aa_colsum[i] > 1 &
			(nt_colsum[nt_position]   > 1 |
			 nt_colsum[nt_position+1] > 1 |
			 nt_colsum[nt_position+2] > 1)) {
			    nt_keep_sites[nt_position] = T
			    nt_keep_sites[nt_position+1] = T
			    nt_keep_sites[nt_position+2] = T
		}
	    }

	    nt_sites <- which(nt_keep_sites)

	} else {

	    nt_sites <- which(nt_colsum > 1)
	}

	R$nts = P$nts[ , nt_sites]
#	R$nt_rast = P$nt_rast[ , nt_sites]
	R$nt_cols = nt_sites
    }

    R
}

#' @keywords internal
create.nj.tree <- function(P) {

    if (class(P) != "pixelgram")
	stop("create.nj.tree ERROR: Please specify pixelgram object")

    message("*** Computing NJ tree from nts using K80, no rate variation ***")

#    S.tree <- ape::bionj(S.dist)
#    S.dist <- ape::dist.dna(ape::as.DNAbin(P$nts))

    P$tre <- ape::bionj(ape::dist.dna(ape::as.DNAbin(P$nts)))

    if (!is.null(P$master_name)) {

	if (!P$master_name %in% c("__consensus__", "__none__")) {
	    P <- try(set.master.name(P))
            P$tre <- ape::root(P$tre, P$master_name)
        }
    }

    P$tre <- ape::ladderize(P$tre)
    P$invert_y = T
    return ( P )
}

# #' @keywords internal
# create.comb.tree <- function(P) {
#
#  if (class(P) != "pixelgram")
#    stop("create.comb.tree ERROR: Please specify pixelgram object")
  # from Will Fischer scripts/branch_color4.pl:
  # generates a bifurcating tree
  # make a "comb" tree (zero-branch-length bifurcations at base)
  # for a list of taxon names
#  @labels      = @_;
#  $newick_tree = pop @labels;
#  while ( my $label = pop @labels ) {
#    newick_tree = paste0("(", label, ":1,", newick_tree, "):1")
#  }
#  newick_tree <- paste0(newick_tree, ';')
# #  newick_tree <- gsub(")", "):0", newick_tree)
#  P$tre <- ape::as.phylo(newick_tree)
#  return (P)
#}

#' @keywords internal
pixelgram.validate <- function(P) {

    if (class(P) != "pixelgram")
	stop("pixelgram.validate ERROR: Please specify pixelgram object")

    message("*** Validating ***")

    if (!P$master_name %in%
          c("__none__", "__consensus__", "__first__", "__last__",
            P$tre$tip.label)) {

        stop(paste0("pixelgram ERROR: Unrecognized master_name: ",
                  P$master_name))
    }

    # TODO (pth): validate nt_master & aa_master
    # if both master name and sequence are set, ensure their equality
    # verify one-to-one mapping of tree tip.labels onto nts

    if (!is.null(P$nts)) {

	if (is.null(P$tre))
	    P <- try(create.nj.tree(P))

        if (length(P$tre$tip.label) != nrow(P$nts))
            cat(paste(
                "pixelgram WARNING: Number of sequences differs between tree and nt alignment: tree has ",
                length(P$tre$tip.label), " leaves vs. ",
                nrow(P$nts), " aligned nts.\n"))

        if (any(!rownames(P$nts) %in% c(P$refseq_name, P$tre$tip.label)))
            stop(paste(
                "pixelgram ERROR: Tree lacks sequences found in nt alignment:",
                paste(rownames(P$nts)[which(!rownames(P$nts) %in% 
		    c(P$refseq_name, P$tre$tip.label))], collapse=" ")))
    }

    # verify one-to-one mapping of tree tip.labels onto aas
    if (!is.null(P$aas)) {

        if (!is.null(P$tre)) {
            if (length(P$tre$tip.label) != nrow(P$aas))
                warning("Number of sequences differs between tree and aa alignment: tree has ",
                    length(P$tre$tip.label), " leaves vs. ",
                    nrow(P$aas), " aligned aas.\n")

            if (any(!rownames(P$aas) %in% c(P$refseq_name, P$tre$tip.label)))
                stop(paste(
                "pixelgram ERROR: Tree lacks sequences found in AA alignment:",
                paste(rownames(P$aas)[which(!rownames(P$aas) %in% 
		    c(P$refseq_name, P$tre$tip.label))], collapse=" ")))
        }
    }

    message("*** Validated OK ***")
    return ( P )
}

#' @keywords internal
pixelgram.reorder <- function(P) {

    if (class(P) != "pixelgram")
	stop("pixelgram.reorder ERROR: Please specify pixelgram object")

    message("*** Reordering ***")

    if (is.null(P$tre))
       return ( P )

    R <- P
    class(R) <- "pixelgram"
    if (!is.null(P$nts)) {

        nt_roworder <- sapply(1:length(P$tre$tip.label), function(i)
            which(rownames(P$nts) == P$tre$tip.label[i]))

        message("*** Reordering nts ***")

        R$nts <- P$nts[nt_roworder, ]
#        R$nt_rast <- P$nt_rast[nt_roworder, ]

        message("*** Reordered nts ok ***")
    }

    if (!is.null(P$aas)) {

        aa_roworder <- sapply(1:length(P$tre$tip.label), function(i)
            which(rownames(P$aas) ==
                  P$tre$tip.label[i]))

        message("*** Reordering aas ***")

        R$aas <- P$aas[aa_roworder, ]
#        R$aa_rast <- P$aa_rast[aa_roworder, ]

        message("*** Reordered aas ok ***")
    }

    R
}

#' @keywords internal
set.master.name <- function(P) {

    if (class(P) != "pixelgram")
	stop("set.master.name ERROR: Please specify pixelgram object")

    if (!is.null(P$master_name)) {

        if (P$master_name %in% c("__first__", "__last__")) {

            # this could all be rewritten elegantly
	    if (!is.null(P$tre)) {
	        if (P$master_name == "__first__") {
		    P$master_name <- P$tre$tip.label[1]
                } else if (P$master_name == "__last__") {
                    P$master_name <- P$tre$tip.label[length(P$tre$tip.label)]
                }
            } else if (!is.null(P$nts)) {
                if (P$master_name == "__first__") {
	            # could/should ensure that refseq has been excised
                    P$master_name <- rownames(P$nts)[1]
                } else if (P$master_name == "__last__") {
                    P$master_name <- rownames(P$nts)[nrow(P$nts)]
                }
            } else if (!is.null(P$aas)) {
                if (P$master_name == "__first__") {

                    P$master_name <- rownames(P$aas)[1]
                } else if (P$master_name == "__last__") {
                    P$master_name <- rownames(P$aas)[nrow(P$aas)]
                }
            }
        }
    }
## set.master.index?
    ## set.refseq?? not needed until rendering in plot()
    return ( P )
}

#' @keywords internal
set.master <- function(P) {

    if (class(P) != "pixelgram")
	stop("set.master ERROR: Please specify pixelgram object")

    if (P$master_name != "__none__") {

        if (P$master_name == "__consensus__") {

            message("*** master is consensus ***")

            if (!is.null(P$nts) & is.null(P$nt_master))
                P$nt_master <- seqinr::consensus(P$nts, "majority")

            if (!is.null(P$aas) & is.null(P$aa_master))
                P$aa_master <- seqinr::consensus(P$aas, "majority")

        } else { # not consensus

	    P <- set.master.name(P)

            if (!is.null(P$nts)) {

                if (!P$master_name %in% rownames(P$nts))
                    stop(paste0("pixelgram ERROR: master '",
                        P$master_name, "' not in nt alignment"))

                if (length(which(rownames(P$nts) == P$master_name)) != 1)
                    stop(paste0("pixelgram ERROR: master '",
                         P$master_name, "' not unique in nt alignment"))

                P$nt_master = P$nts[P$master_name, ]
	    }

            if (!is.null(P$aas)) {

                if (!P$master_name %in% rownames(P$aas))
                    stop(paste0("pixelgram ERROR: did not find master '",
                        P$master_name, "' in aa alignment"))

                if (length(which(rownames(P$aas) == P$master_name)) != 1)
                    stop(paste0("pixelgram ERROR: master '",
                        P$master_name, "' not unique in aa alignment"))

                P$aa_master = P$aas[P$master_name, ]
	    }
	    message(paste0("*** master is ", P$master_name))
        }
    }
    return ( P )
}
#' @keywords internal
pixelgram.xform <- function(P) {

    if (class(P) != "pixelgram")
	stop("pixelgram.xform ERROR: Please specify pixelgram object")

    if (!P$xform_type %in% c(0:3))
	stop("pixelgram.xform ERROR: Please indicate transformation type.")

    if (P$master_name == "__none__") {
	P$aa_rast = P$aas
	P$nt_rast = P$nts
	return ( P )
    }
    if (!is.null(P$aas) & !is.null(P$aa_master)) {

        message("*** Xforming AAs ***")

	P$aa_rast = P$aas

        if (length(P$aa_master) != ncol(P$aas))
            stop(paste("pixelgram ERROR: aa master mismatches alignment length:",
		    length(P$aa_master), '!=', ncol(P$aas)))

        # if master is in the alignment, need to restore it after xforming
        master_index <- ifelse(P$master_name == "__consensus__", 0,
            which(rownames(P$aas)==P$master_name))

#	P$aa_master <- P$aas[master_row, ]

        if (P$xform_type == 1) {

            # =:matches(white), %:indels(black), !:mismatches(red)
            # $s_first eq $s_curr ? '=' :
            #     (($s_first eq '-' or $s_curr eq '-') ? '%' : '!')

            for (i in 1:ncol(P$aas))
                # all sites that match the master are '='
                # gain or loss of "-" should become '%'
                # point mutations are "!"
                P$aa_rast[, i] <- gsub('[a-y]', '!',       # 3: mark intact codons
                    replace(replace(P$aas[ ,i],  # 1: change matches to master
			    which(P$aas[ ,i] == P$aa_master[i]), '='),
			which(xor(P$aa_master[i] == '-', P$aas[, i] =='-')), '%')) # 2: insertions
# to do: consider other colors for inviable codons
        } else if (P$xform_type == 2) {

            # @:matches(grey) -> = maybe . instead
            # =:deletion(white) -> %
            # original base:mismatches(usual color)
            # $s_first eq $s_curr ? '@' : ( $s_curr eq '-' ? '=' : $s_curr)

            for (i in 1:ncol(P$aas))
                P$aa_rast[, i] <- gsub('-', '%', replace(P$aas[ ,i],
                    which(P$aas[ ,i] == P$aa_master[i]), '.'))

        } else if (P$xform_type == 3) {

#[^RKH]->[RKH]
#[^ED]->[ED]
#[^O]->O

#        message("    Please be patient while computing this transformation.")

	misc.aa <- c("!", "#", "$", "z", "x", "+", "_", "-", "%", "=", "o")
	pos.aas <- c("k", "r", "h")
	neg.aas <- c("d", "e")

        for (i in 1:ncol(P$aas)) {

            x <- replace(P$aas[, i],
		which(P$aas[, i] == P$aa_master[i]), '=')

	    if (!P$aa_master[i] %in% pos.aas)
		x <- replace(x, which(x %in% pos.aas), "+")

	    if (!P$aa_master[i] %in% neg.aas)
		x <- replace(x, which(x %in% neg.aas), "_")

	    if (P$aa_master[i] != "-")
		x <- replace(x, which(x == "-"), "%")

	    P$aa_rast[, i] <- replace(x, which(!x %in% misc.aa), "^")

	}

# this was inefficient:
#                  P$aa_rast[, i] <- sapply(1:nrow(P$aas), function(j)
#  		    ifelse(P$aas[j,i] == P$aa_master[i], '=', # match t/f
#  			ifelse(P$aas[j,i] == 'o', 'o', # gain of png site
#  			    ifelse(!P$aa_master[i] %in% c('r','k','h') &
#                               P$aas[j,i] %in% c('r','k','h'), '+',
#  				    ifelse(!P$aas[i] %in% c('d','e') &
#                                       P$aas[j,i] %in% c('d','e'), '_',
#   				            ifelse(P$aas[j,i] == '-' |
#  					    P$aa_master[i]=='-', '%','^'))))))
        }

	if (!P$xform_master & master_index) # if master_index == 0 (consensus) don't restore it
	    P$aa_rast[master_index, ] <- P$aa_master
    }

    if (!is.null(P$nts) & !is.null(P$nt_master)) {

        message("*** Xforming NTs***")

	P$nt_rast = P$nts

        if (length(P$nt_master) != ncol(P$nts))
            stop(paste("pixelgram ERROR: nt master mismatches alignment:",
		length(P$nt_master), '!=', ncol(P$nts)))

        # if master is in the alignment, need to restore it after xforming
        master_index <- ifelse(P$master_name == "__consensus__", 0,
            which(rownames(P$nts) == P$master_name))

        if (P$xform_type == 1) {

# pixel code does this:
# =:matches(white), %:indels(black), !:mismatches(red)
# $s_first eq $s_curr ? '=' : (($s_first eq '-' or $s_curr eq '-') ? '%' : '!')

            # all sites that match the master are '='
            for (i in 1:ncol(P$nt_rast))
                P$nt_rast[, i] <- gsub('-', '%', gsub('[a-z]', '!',
                    replace(P$nts[ ,i],
                        which(P$nts[ ,i] == P$nt_master[i]),'=')))

# annotate hypermutations
	    for (i in which(P$nt_master[1:(length(P$nt_master)-2)] == 'g'))
                P$nt_rast[which(P$nts[, i  ] == 'a'
			     & (P$nts[, i+1] == 'a' |
				P$nts[, i+1] == 'g')
			     & (P$nts[, i+2] == 'a' |
			        P$nts[, i+2] == 'g' |
				P$nts[, i+1] == 't') ), i ] = '^'
        } else if (P$xform_type > 1) {

# pixel code says:
# @:matches(grey), =:deletion(white), original base:mismatches(usual color)
# $s_first eq $s_curr ? '@' : ( $s_curr eq '-' ? '=' : $s_curr)

            for (i in 1:ncol(P$nt_rast))
                P$nt_rast[, i] <- gsub('-', '%',
                    replace(P$nts[ ,i],
                        which(P$nts[ ,i] == P$nt_master[i]), '='))
	    }

            if (!P$xform_master & master_index) # non-zero row
                P$nt_rast[master_index, ] <- P$nt_master
	}

     message("*** Xformed ***")

    return ( P )
}

#' @export
print.pixelgram <- function(x, ...) summary.pixelgram(x, ...)

#' @export
summary.pixelgram <- function(object, ...) {

    if (class(object) != "pixelgram")
	stop("pixelgram.summary ERROR: Please specify pixelgram object")

    message(paste("Tree has", length(object$tre$tip.label), "leaves."))

    if (!is.null(object$nts))
        message(paste0("Aligned nts: ", dim(object$nts)[1],
                     " sequences, length ", dim(object$nts)[2]))

    if (!is.null(object$aas))
        message(paste0("Aligned aas: ", dim(object$aas)[1],
                     " sequences, length ", dim(object$aas)[2]))

    if (!is.null(object$master_name))
        message(paste0("Using '", object$master_name, "' for master sequence."))
}
