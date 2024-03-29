% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/linreg.R
\docType{class}
\name{linreg-class}
\alias{linreg-class}
\alias{linreg}
\title{this class is calculating Regressions coefficients, fitted values, residuals, 
degrees of freedom, residual variance, variance of the regression coefficients,
t-values for each coefficient and p-value for each coefficient.}
\arguments{
\item{formula}{this is the formula that need to be passed into contructor of the class, it will contain variable names from data frame.}

\item{data}{this will be the data frame that we are working on.}
}
\value{
this will return class object, which will then be used to access multiple function of the class.
}
\description{
this class is calculating Regressions coefficients, fitted values, residuals, 
degrees of freedom, residual variance, variance of the regression coefficients,
t-values for each coefficient and p-value for each coefficient.
}

\examples{
\dontrun{
data("iris")
mat_obj<-linreg$new()
mat_obj<-mat_obj$linreg(formula=Petal.Length~Species,data = iris)
mat_obj$print()
mat_obj$plot()
mat_obj$resid()
mat_obj$summary()
}
}
