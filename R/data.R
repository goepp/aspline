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
#' @source \code{fda::"MontrealTemp"}
"montreal"
#' Fossil data
#'
#' A dataset with 106 observations on fossil shells from the \code{SemiPar} package (\url{https://CRAN.R-project.org/package=SemiPar}).
#'
#' @format A data frame with 106 observations and 2 variables:
#' \describe{
#'   \item{age}{The age of fossils, in millions of years}
#'   \item{strontium.ratio}{Ratio of strontium isotopes}
#'   ...
#' }
#' @source Bralower, T.T, Fullagar, P.D., Paull, C.K, Dwyer, G.S. and Leckie, R.M. (1997). Mid-cretaceous strontium-isotope stratigraphy of deep-sea sections. \emph{Geological Society of America Bulletin}, 109, 1421-1442.
#' @references Ruppert, D., Wand, M.P. and Carroll, R.J. (2003). \emph{Semiparametric Regression}, Cambridge University Press.
"fossil"
#' Bladder Cancer aCGH profile data
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
#' of Transcription in Cancer', \emph{Nature Genetics} 38(12), 1386-1396.
"bladder"
#' Yearly number of coal mine disasters in Britain
#'
#' A data of 112 observations registering the yearly number of coal mine disasters in Britain from 1851 to 1962.
#' The data comes from \cite{Diggle et al. (1988)} and has been used for spline regression by \cite{Eilers et al. (1996)}.
#'
#' @format A data frame with 112 observations and 2 variables:
#' \describe{
#'   \item{year}{year}
#'   \item{n}{number of coal mine disasters}
#' }
#' @source Diggle, P. and Marron, J. S. (1988). `Equivalence of Smoothing Parameter Selectors
#' in Density and Intensity Estimation', \emph{Journal of the American Statistical Association} 83(403), 793-800.
#' @references Eilers, P. H. C. and Marx, B. D. (1996). `Flexible Smoothing with B-splines and Penalties', \emph{Statistical Science} 11(2), 89-102.
"coal"
#' Titanium heat data
#'
#' A data set of 49 samples expressing the thermal property of titanium
#'
#' @format 49 observations and two variables:
#' \describe{
#' \item{x}{temperature}
#' \item{y}{physical property}
#' }
#' @source
#' \itemize{
#' \item de Boor, C., and Rice, J. R. (1986), Least-squares cubic spline approximation. II: variable knots. \emph{Report CSD TR 21, Purdue U., Lafayette, IN.}
#' \item Dierckx, P. (1993), \emph{Curve and Surface Fitting with Splines}, Springer.
#' \item Jupp, D. L. B. (1975), \emph{Approximation to data by splines with free knots}, SIAM Journal on Numerical Analysis, 15: 328-343.
#' }
#'
"titanium"
#' Nuclear Magnetic Resonance data
#'
#' A signal of nuclear magnetic resonance.
#'
#' @format Data farme of 1024 rows and two columns: the index \code{x} and the signal \code{y}.
#' @source
#' \itemize{
#' \item Data from \url{https://web.stanford.edu/~hastie/ElemStatLearn/datasets/nmr1.csv}.
#' \item See also The Elements of Statisical Learning (2001, 2nd Ed.), \emph{Hastie, T., Friedman, J., and Tibshirani, R.J}, p. 176}.
#'
"nmr"
#' Lidar data
#'
#' Data from a light detection and ranging (LIDAR) experiment
#'
#' @format
#' \describe{
#'    \item{range}{distance travelled before the light is reflected back to its source}
#'    \item{logratio}{logarithm of the ratio of received light from two laser sources}}
#' @source
#' \itemize{
#' \item Sigrist, M. (Ed.) (1994). Air Monitoring by Spectroscopic Techniques (Chemical Analysis Series,vol. 197). New York: Wiley
#' \item The R package \url{https://CRAN.R-project.org/package=SemiPar}}
#' @references Ruppert, D., Wand, M.P. and Carroll, R.J. (2003). \emph{Semiparametric Regression}, Cambridge University Press.
"lidar"
