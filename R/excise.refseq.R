#' @keywords internal
excise.refseq <- function(refseq_row, aln) {

    # return alignment wihthout aas_aln[refseq_row, ] has one less row
    message(paste0("Excising refseq (row ", refseq_row, ") from alignment."))
    aln[-refseq_row, ]

  # this is invoked by set.refseq.R
}
