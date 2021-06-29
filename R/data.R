#' Ensemble forecasts of temperature and associated observation
#'
#' This is a data set containing the forecasts of five ensemble weather
#' prediction models for two-meter temperature from March, 2019 to March, 2021.
#'
#' The five ensemble models are named CWAO (20 perturbed members, from ECCC),
#' DEMS (11 perturbed members, from NCMRWF), ECMF (50 members, from ECMWF),
#' EGRR (17 perturbed members, from UKMO) and RKSL (24 perturbed members, from
#' KMA).
#'
#' The forecasts are the ones at the nearest grid point to Toulouse-Blagnac
#' station (France) in the TIGGE dataset. The observation is the two-meter
#' height temperature measured at this same station, at 06UTC. The forecast
#' initial time is 00UTC, with a 30 hour lead-time.
#'
#' @format A named list with five entries, each containing a \code{data.table}
#' with 731 rows and a varying number of columns (depending on the number of
#' members).
#' \describe{
#' \item{date_run}{initial condition date (with format YYYY-MM-DD)}
#' \item{latitude}{latitude of the forecast location (in degrees)}
#' \item{longitude}{longitude of the forecast location (in degrees)}
#' \item{1 ... M}{member forecast (M members) in Kelvins}
#' \item{T}{the measured two-meter air temperature in Kelvins}
#' }
#' @source  \url{https://apps.ecmwf.int/datasets/data/tigge/levtype=sfc/type=pf/}
#'
#' \url{https://donneespubliques.meteofrance.fr/?fond=produit&id_produit=91&id_rubrique=32}
"ensembles"

#' Post-processed ensemble forecasts of temperature and associated observation
#'
#' This is a dataset containing the post-processed forecasts of five ensemble
#' weather prediction models for two-meter temperature from March, 2019 to
#' March, 2021.
#'
#' Each ensemble has been post-processed with the non homogeneous regression
#' technique, described in Gneiting et al.(2005). In a nutshell the true
#' distribution is supposed to be gaussian, with mean and standard deviation
#' being a linear function of the ensemble mean and standard deviation
#' (respectively). The intercept and slope of each regression is determined by
#' minimizing the CRPS over a 60-day sliding window. The forecast in the
#' data set is a sample of 30 values from this gaussian distribution.
#'
#' The five ensemble models are named based on the raw ensemble: CWAO
#' (from ECCC), DEMS (from NCMRWF), ECMF (from ECMWF), EGRR (from UKMO) and RKSL
#' (from KMA).
#'
#' The raw forecasts are the ones at the nearest grid point to Toulouse-Blagnac
#' station (France) in the TIGGE data set. The observation is the two-meter
#' height temperature measured at this same station, at 06UTC. The forecast
#' initial time is 00UTC, with a 30 hour lead-time.
#'
#' @format A named list with five entries, each containing a \code{data.table}
#' with 731 rows.
#' \describe{
#' \item{date_run}{initial condition date (with format YYYY-MM-DD)}
#' \item{latitude}{latitude of the forecast location (in degrees)}
#' \item{longitude}{longitude of the forecast location (in degrees)}
#' \item{1 ... 30}{forecast in Kelvins sampled from the gaussian distribution.
#' The forecasts are sorted in increasing order.}
#' \item{T}{the measured two-meter air temperature in Kelvins}
#' }
#' @source  \url{https://apps.ecmwf.int/datasets/data/tigge/levtype=sfc/type=pf/}
#'
#' \url{https://donneespubliques.meteofrance.fr/?fond=produit&id_produit=91&id_rubrique=32}
#' @references
#' Gneiting, Tilmann, et al. "Calibrated probabilistic forecasting using
#' ensemble model output statistics and minimum CRPS estimation."
#' \emph{Monthly Weather Review} 133.5 (2005): 1098-1118.
#' doi:https://doi.org/10.1175/MWR2904.1
"ppensembles"
