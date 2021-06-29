## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  options(rmarkdown.html_vignette.check_title = FALSE)
)

## -----------------------------------------------------------------------------
set.seed(42)
N <- 20
K <- 21
ranks <- sample(1:K, size = K * N, replace = TRUE)
shst <- hst <- hist(ranks, breaks = 0:K + 0.5)
abline(h = N, lty = 2)
shst$counts <- sort(shst$counts)
plot(shst, col ="lightgrey")
abline(h = N, lty = 2)
cat("Random counts' chisquare statistics =", chisq.test(hst$counts)$statistic, "\n")
cat("Sorted counts' chisquare statistics =", chisq.test(shst$counts)$statistic, "\n")

## -----------------------------------------------------------------------------
library(flatness)
ensembles <- ensembles
is(ensembles)
names(ensembles)

## -----------------------------------------------------------------------------
library(data.table)
ensembles$CWAO
range(ensembles$CWAO$date_run)
sapply(ensembles, function(x) NCOL(x) - 4)

## -----------------------------------------------------------------------------
ppensembles <- ppensembles
print(names(ppensembles))
print(str(ppensembles$ECMF))
sapply(ppensembles, function(x) NCOL(x) - 4)

## -----------------------------------------------------------------------------
iobs <- which(names(ensembles$CWAO) == "T")
ifcst <- 4:(iobs - 1)
M <- length(ifcst)
K <- M + 1
print(paste("Number of members:", M))
fcst <- data.matrix(ensembles$CWAO[, ifcst, with = FALSE])
obs <- unlist(ensembles$CWAO[, iobs, with = FALSE])
rkh <- rkhist(fcst = fcst, obs = obs)
class(rkh)
str(rkh)

## -----------------------------------------------------------------------------
plot(rkh)
plot(rkh, what = "percents")
rownames(rkh) <- "CWAO"
print(rkh)
print(rkh, what = "percents", digits = 1L)

## -----------------------------------------------------------------------------
linear <- deviate_linear(k = K)
plot(linear)
U <- deviate_U(k = K)
plot(U)
V <- deviate_V(k = K)
plot(V)
ends <- deviate_ends(k = K)
plot(ends)
wave <- deviate_wave(k = K)
plot(wave)
print(is(wave))

## -----------------------------------------------------------------------------
deviates <- get_deviates(k = K, shapes = c("linear", "U", "V", "ends", "wave"))
print(str(deviates))
plot(deviates)

## -----------------------------------------------------------------------------
check <- is_JP_ready(x = deviates, tol = 0.001)
print(str(check))

## -----------------------------------------------------------------------------
deviates_ok <- make_JP_ready(x = deviates)
check <- is_JP_ready(x = deviates_ok)
plot(deviates_ok)

## -----------------------------------------------------------------------------
deviates <- get_deviates(k = K)
print(row.names(deviates))

## -----------------------------------------------------------------------------
JP <- JP_test(rkhists = rkh, deviates = deviates)
print(str(JP))

## -----------------------------------------------------------------------------
print(JP)
print(JP, what = "pvalues")
library(lattice)
levelplot(JP)
levelplot(JP, what = "projections", main = "projections")

## -----------------------------------------------------------------------------
ppensembles <- ppensembles
fcsts <- lapply(ppensembles, function(x) x[, 4:33, with =FALSE])
obss <- lapply(ppensembles, function(x) t(x[, 34, with =FALSE]))
rkhists <- rkhist(fcsts, obss)
print(rkhists)
plot(rkhists, what = "percents")
plot(rkhists)
deviates <- get_deviates(k = NCOL(rkhists))
ppJP <- JP_test(rkhists, deviates)
levelplot(ppJP, main = "pvalues")
levelplot(ppJP, what = "projections", main = "projections")

## -----------------------------------------------------------------------------
indices <- flatness_indices(rkhists)
print(indices)

