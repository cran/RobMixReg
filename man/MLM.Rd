% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/outer2.R
\name{MLM}
\alias{MLM}
\title{The main function of mining the latent relationship among variables.}
\usage{
MLM(ml.method = "rlr", rmr.method = "cat",
  b.formulaList = list(formula(y ~ x), formula(y ~ 1)), formula = y ~
  x, nit = 1, nc = 2, x = NULL, y = NULL, max_iter = 50,
  tRatio = 0.05)
}
\arguments{
\item{ml.method}{The option to select the four methods in vignette.}

\item{rmr.method}{The option to select the robust mixture regression method.}

\item{b.formulaList}{The case b require the user provide the formula list. This enable the flexible mixture regression.}

\item{formula}{The linear relationship between two variables.}

\item{nit}{Number of iterations for CTLE, mixbi, mixLp.}

\item{nc}{Number of mixture components.}

\item{x}{The matrix x of the high dimension situation.}

\item{y}{The external outcome variable.}

\item{max_iter}{Maximum iteration for TLE method.}

\item{tRatio}{The ratio of the outliers in the TLE robust mixture regression method.}
}
\value{
Main result object.
}
\description{
The main function of mining the latent relationship among variables.
}
