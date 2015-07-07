#' Set the output file name and format for saving a rendered pixgram.
#'
#' @param pg Pixgram object.
#' @param outfile_name If provided, the output will be written to this file name.
#' @param outfile_format If outfile_name is null, this is ignored and the output is written directly with plot(); otherwise, this must be one of "png", "svg", "eps", "pdf", in which case the device-specific driver is called before plot() and closed afterwards.
#'
#' @export
set.output.file <- function(pg, outfile_name, outfile_format) {

    if (class(pg) != "pixgram")
	stop("set.output.file ERROR: Please specify pixgram object")

    if (is.null(outfile_name)) {
        warning("Null file name")
	return ( pg )
    }

    if (!outfile_format %in% c("png", "svg", "eps", "pdf")) {
	warning("Invalid file format")
	return ( pg )
    }    

    pg$outfile_format = outfile_format
    pg$outfile_name = outfile_name

    switch ( pg$outfile_format,
	    png = png(filename=pg$outfile_name),
	    svg = svg(filename=pg$outfile_name),
	    eps = postscript(file=pg$outfile_name),
	    pdf = pdf(file=pg$outfile_name, useDingbats=F) )

    par(mar=c(0,0,0,0), oma=c(0,0,0,0), lend=1, ljoin=1)
    return ( pg )
}
