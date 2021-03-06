% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mixL.short.R
\name{mixLp_one}
\alias{mixLp_one}
\title{mixLp_one : mixLp_one estimates the mixture regression parameters robustly using Laplace distribution based on one initial value.}
\usage{
mixLp_one(formula, data, nc = 2)
}
\arguments{
\item{formula}{A symbolic description of the model to be fit.}

\item{data}{A data frame containing the predictor and response variables, where the last column is the response varible.}

\item{nc}{Number of mixture components.}
}
\value{
Estimated coefficients of all components.
}
\description{
Robust mixture regression assuming that the error terms follow a Laplace distribution.
}
