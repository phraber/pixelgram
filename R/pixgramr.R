#' pixgramr: Display a pixel plot and phylogram of aligned sequences.
#'
#' The two principal methods are aconstructor, which usies the S3 object model and  a plotting function.  For details about how to use each, run \code{help(pixgramr::pixgram)} and \code{help(pixgramr::plot.pixgram)}.
#'
#' For examples run \code{vignette("pixgramr")}.
#'
#' @docType package
#' @name pixgramr
#' @seealso \code{\link{pixgram}} \code{\link{plot}}
NULL

#' @export
pixgram <- function(...) UseMethod("pixgram")
