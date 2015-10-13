#' @keywords internal
get.refseq.row <- function(aln, refseq_name, ignore_case=F) {

    if (is.null(refseq_name))
	return(NULL)

    name.hits <- grepl(refseq_name, rownames(aln), ignore.case=ignore_case)
    n.matches <- length(which(name.hits))

    if (n.matches == 1) {

	return(which(name.hits))

    } else if (n.matches == 0) {

	if (ignore_case) { # prevent infinite regress
	    return(NULL)
	} else {
	    get.refseq.row(aln, refseq_name, ignore_case=T)
	}

    } else {
	stop(paste0("pixelgram::get.refseq.row() ERROR: refseq_name '",
	    refseq_name, "' found more than once in alignment.\n",
	    "To fix this, please delete the extra entries.\n"))
    }
}

