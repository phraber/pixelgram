#' @keywords internal
create.refseq.lut <- function(P, ref) {

    if (class(P) != "pixelgram")
        stop("ERROR: Please pass a pixelgram object to create.refseq.lut()")

    if (is.null(P$refseq_lut) & !is.null(ref)) {

	aln <- c(1:length(ref))

	l <- rep(NA, length(ref))
	r <- rep(NA, length(ref))

	L <- 0
	R <- 1

    # there must be an efficient alternative!
	for (i in 1:length(ref)) {

	    if (ref[i] != "-") L = L + 1

		l[i] = L
		r[i] = R

		if (ref[i] != "-") R = R + 1
	}

	P$refseq_lut <- data.frame(aln, l, r, ref)
    } # don't change the lut if it is already set

    return ( P )
}
