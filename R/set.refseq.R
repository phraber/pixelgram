#' @keywords internal
set.refseq <- function(P, is.aa=NULL, use.grepl=F) {

    if (class(P) != "pixgram")
        stop("ERROR: Please pass a pixgram object to set.refseq()")

    if (is.aa==F) {
        aln <- P$nts
    } else if (is.aa==T) {
        aln <- P$aas
    } else {
	return ( P )
    }

    if (is.null(P$refseq_row) & !is.null(P$refseq_name)) {
    # set refseq_row by matching refseq_name

	if (use.grepl) { # use only if an exact match is not found first-time through; see below
	    if (length(which(grepl(P$refseq_name, rownames(aln)))) == 1) {
	        P$refseq_name = rownames(aln)[which(grepl(P$refseq_name, 
		    rownames(aln)))]

		P$refseq_row <- get.refseq.row(aln, P$refseq_name)
	    }
        } else {

	    if (P$refseq_name %in% rownames(aln))
		P$refseq_row <- get.refseq.row(aln, P$refseq_name)
	}
    }

    ### create lut and excise refseq from alignment if it's in there
    if (!is.null(P$refseq_row)) {

	P <- create.refseq.lut(P, aln[P$refseq_row, ])

        ### HERE WE EXCISE REFSEQ from ALIGNMENT:
	if (!is.null(aln) & P$excise_refseq & !is.null(P$refseq_name) & 
	    P$refseq_name == rownames(aln)[P$refseq_row] & !use.grepl) # only if names exactly match
	        aln <- excise.refseq(P$refseq_row, aln)

	if (!is.null(is.aa))
	    if (!is.aa)
		if (!is.null(P$nts))
		    P$nts = aln

	if (!is.null(is.aa))
	    if (is.aa)
		if (!is.null(P$aas))
		    P$aas = aln

    } else if (is.null(P$refseq_row)) {


# try again with regexps, e.g. HXB2 in full name
	if (!use.grepl)
	    set.refseq(P, is.aa, use.grepl=T)

        # no reference sequence - just use alignment master for numbering
#        if (!is.null(P$master_index)) {
#	    P$refseq_lut <- create.refseq.lut(P, P$aas_aln[P$master_index, ])
#	    P$refseq_name <- rownames(P$aas_aln)[P$master_index]
#	    P$refseq_row <- P$master_index
#        }
    }

    return ( P )
}

