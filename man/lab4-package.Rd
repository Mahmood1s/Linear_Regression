\name{lab4-package}
\alias{lab4-package}
\alias{lab4}
\docType{package}
\title{
\packageTitle{lab4}
}
\description{
\packageDescription{lab4}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{lab4}
\packageIndices{lab4}
~~ this package is calculating the linear regression on the dataset provided and the fomula passed to the class. package is using RC class implementation for linear Regression. there are various function which are overloaded like print, summary, plot, coef, pred and resid.  ~~
}
\author{
\packageAuthor{lab4}

Maintainer: \packageMaintainer{lab4}
}
\examples{
\dontrun{
mat_obj<-linreg$new(Petal.Length~Sepal.Width+Sepal.Length,data = iris)
mat_obj$print()
mat_obj$plot()
mat_obj$pred()
mat_obj$resid()
mat_obj$summary()
}
}
