#' Example dataset based on cross-validated prediction of coronary disease
#' in the Cleveland Heart Study 
#'
#' @format A data frame with 768 rows and three variables
#' \describe{
#'   \item{prior.p}{prior probabilities of case status}
#'   \item{posterior.p}{posterior probabilities of case status}
#'   \item{y}{case-control status}
#' }
#' @source \url{http://www.homepages.ed.ac.uk/pmckeigu/preprints/classify/demoplotw.html}
"cleveland"
#'
#' Example dataset based on cross-validated prediction of diabetes
#' in Pima Native Americans
#'
#' @format A data frame with 297 rows and three variables
#' \describe{
#'   \item{prior.p}{prior probabilities of case status}
#'   \item{posterior.p}{posterior probabilities of case status}
#'   \item{y}{case-control status}
#' }
#' @source \url{http://www.homepages.ed.ac.uk/pmckeigu/preprints/classify/demoplotw.html}
"pima"
#'
#' Example dataset based on cross-validated prediction of colorectal cancer 
#' from fecal immunochemical test (FIT) only in Michigan.  As most controls and
#' some cases have have zero values in the FIT test, to fit densities to the sampled
#' values of weight of evidence in controls and cases it is necessary to specify spike-slab
#' mixtures. 
#'
#' @format A data frame with 242 rows and three variables
#' \describe{
#'   \item{prior.p}{prior probabilities of case status}
#'   \item{posterior.p}{posterior probabilities of case status}
#'   \item{y}{case-control status}
#' }
#' @source \url{http://www.homepages.ed.ac.uk/pmckeigu/preprints/classify/demoplotw.html}
"fitonly"
#'
#' Example dataset based on cross-validated prediction of colorectal cancer
#' from fecal immunochemical test (FIT) plus microbiome profile in Michigan
#'
#' @format A data frame with 242 rows and three variables
#' \describe{
#'   \item{prior.p}{prior probabilities of case status}
#'   \item{posterior.p}{posterior probabilities of case status}
#'   \item{y}{case-control status}
#' }
#' @source \url{http://www.homepages.ed.ac.uk/pmckeigu/preprints/classify/demoplotw.html}
"fitplusmicrobiome"
#'
 


