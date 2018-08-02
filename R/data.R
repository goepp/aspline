#' Testing Crash Helmets
#'
#' A dataset containing the acceleration and time after impact of helmets
#' from a simulated motorcycle accident.
#'
#' @format A data frame with 132 rows and 2 variables:
#' \describe{
#'   \item{x}{Time after impact, in milliseconds}
#'   \item{y}{Head acceleration, in units of \eqn{g}}
#'   ...
#' }
#' @source Dataset number \eqn{338} of \emph{Hand, D. et al. (1993) A Handbook of Small Datasets}.
"helmet"
#' Montreal Temperature Data
#'
#' A dataset containing the tempature in Montreal for two years
#'
#' @format A data frame with 730 rows and 2 variables:
#' \describe{
#'   \item{day}{The day of the year from January 1, 1961, to December 31, 1962}
#'   \item{temp}{Temperature in Celsius}
#'   ...
#' }
#' @source \url{http://people.stat.sfu.ca/~cschwarz/Consulting/Trinity/Phase2/TrinityWorkshop/Workshop-material-Simon/Intro_to_splines/Exercises/}
"montreal"
#' Fossil data
#'
#' A dataset with 106 observations on fossil shells from the \code{SemiPar} package (\url{https://CRAN.R-project.org/package=SemiPar}).
#'
#' @format A data frame with 106 observations and 2 variables:
#' \describe{
#'   \item{age}{The age of fossils, in millions of years}
#'   \item{strontium}{Ratio of strontium isotopes}
#'   ...
#' }
#' @source Bralower, T.T, Fullagar, P.D., Paull, C.K, Dwyer, G.S. and Leckie, R.M. (1997). Mid-cretaceous strontium-isotope stratigraphy of deep-sea sections. \italic{Geological Society of America Bulletin}, 109, 1421-1442.
#' @references Ruppert, D., Wand, M.P. and Carroll, R.J. (2003). \italic{Semiparametric Regression}, Cambridge University Press.
#' \url{http://stat.tamu.edu/~carroll/semiregbook/}
"fossil"
#' Bladder Cancer aCGH data
#'
#' A dataset of 500 observations corresponding to 500 probes of the aCGH profile of a bladder cancer patient.
#' The original data are provided by \cite{Stransky et al. (2006)}. This dataset consists of probes 1 through 500 of individual 1.
#'
#' @format A data frame with 500 observations and 2 variables:
#' \describe{
#'   \item{x}{probe number}
#'   \item{y}{aCGH profile value}
#' }
#' @source Stransky, N., Vallot, C., Reyal, F., Bernard-Pierrot, I., de Medina, S. G. D., Segraves,
#' R., de Rycke, Y., Elvin, P., Cassidy, A., Spraggon, C., Graham, A., Southgate, J.,
#' Asselain, B., Allory, Y., Abbou, C. C., Albertson, D. G., Thiery, J. P., Chopin, D. K.,
#' Pinkel, D. and Radvanyi, F. (2006). Regional Copy Number Independent Deregulation
#' of Transcription in Cancer', \italic{Nature Genetics} 38(12), 1386-1396.
"bladder"
