#' Example datasets
#'
#' The \pkg{wevid} package comes with the following dataset:
#' \itemize{
#' \item \code{cleveland} is based on cross-validated prediction of coronary
#'       disease in the Cleveland Heart Study (297 observations).
#' }
#'
#' @name wevid.datasets
#' @docType data
#' @usage NULL
#' @format Each dataset consists of a data frame with the following variables:
#' \itemize{
#'   \item prior.p: Prior probabilities of case status.
#'   \item posterior.p: Posterior probabilities of case status.
#'   \item y: Case-control status.
#' }
#' @source \url{http://www.homepages.ed.ac.uk/pmckeigu/preprints/classify/wevidtutorial.html}
#' @keywords datasets
"cleveland"

#' \itemize{
#' \item \code{pima} is based on cross-validated prediction of diabetes
#'       in Pima Native Americans (768 observations).
#' }
#'
#' @usage NULL
#' @rdname wevid.datasets
"pima"

#' \itemize{
#' \item \code{fitonly} is based on cross-validated prediction of colorectal
#'       cancer from fecal immunochemical test (FIT) only in Michigan (242
#'       observations). As most controls and some cases have have zero values
#'       in the FIT test, to fit densities to the sampled values of weight of
#'       evidence in controls and cases it is necessary to specify spike-slab
#'       mixtures.
#' }
#'
#' @usage NULL
#' @rdname wevid.datasets
"fitonly"
