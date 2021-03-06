% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/outer.R
\name{rmr}
\alias{rmr}
\title{The main function of Robust Mixture Regression using five methods.}
\usage{
rmr(lr.method = "flexmix", formula = NULL, data = NULL, nc = 2,
  nit = 20, tRatio = 0.05, MaxIt = 200)
}
\arguments{
\item{lr.method}{A robust mixture regression method to be used. Should be one of "flexmix", "TLE", "CTLERob", "mixbi","mixLp".}

\item{formula}{A symbolic description of the model to be fit.}

\item{data}{A data frame containing the predictor and response variables, where the last column is the response varible.}

\item{nc}{Number of mixture components.}

\item{nit}{Number of iterations for CTLE, mixbi, mixLp.}

\item{tRatio}{Trimming proportion for TLE method.}

\item{MaxIt}{Maximum iteration for TLE method.}
}
\value{
An S4 object about the regression result.
}
\description{
The main function of Robust Mixture Regression using five methods.
}
\examples{
library(RobMixReg)
#library(robust)
library(flexmix)
library(robustbase)
library(MASS)
library(gtools)
# gaussData
x=(gaussData$x);y=as.numeric(gaussData$y);
formula01=as.formula("y~x")
example_data01=data.frame(x,y)
res_rmr = rmr(lr.method='flexmix', formula=formula01, data=example_data01)
res_rmr = rmr(lr.method='CTLERob', formula=formula01, data=example_data01)

}
