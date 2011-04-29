\encoding{latin1}
\name{plot.valid}
\alias{plot.valid}

\title{Validation Plot}

\description{
Function to plot validation statistics, such as MSEP, RMSEP or \eqn{R^2}, as a
function of the number of components.
}

\usage{
\method{plot}{valid}(x, criterion = c("MSEP", "RMSEP", "R2"),
     pred.method = "all",
     xlab = "number of components", ylab = NULL,
     cTicks = NULL, layout = NULL, \dots)
}

\arguments{
  \item{x}{an \code{valid} object.}
  \item{criterion}{character string. What type of validation criterion to plot for 
    \code{pls} or \code{spls}. One of \code{"MSEP"}, \code{"RMSEP"} or \code{"R2"}. 
	See \code{\link{valid}}.}
  \item{pred.method}{prediction method applied in \code{valid} for \code{plsda} or \code{splsda}. 
	See \code{\link{valid}}.}
  \item{xlab, ylab}{titles for \eqn{x} and \eqn{y} axes.  Typically character strings, 
    but can be expressions (e.g., \code{expression(R^2)}).}  
  \item{cTicks}{integer vector. Axis tickmark locations for the used number of components.
    Default is \code{1:ncomp} (see \code{\link{valid}}).}
  \item{layout}{numeric vector of length two giving the number of rows and columns in a 
    multi panel display. If not specified, \code{plot.valid} tries to be intelligent.}
  \item{\dots}{Further arguments sent to \code{\link{xyplot}} function.}
}

\details{
\code{plot.valid} creates one plot for each response variable in the
model, laid out in a multi panel display.  It uses \code{\link{xyplot}} for
performing the actual plotting.
}

\author{Ignacio González.}

\seealso{\code{\link{pls}}, \code{\link{spls}}, \code{\link{predict}}.}

\examples{
require(lattice)

## validation for objects of class 'pls' or 'spls'
\dontrun{
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic

liver.val <- valid(X, Y, ncomp = 5, mode = "regression", 
                   method = "pls", validation = "Mfold")
				   
plot(liver.val, criterion = "R2", type = "l", layout = c(2, 2))
}

## validation for objects of class 'plsda' or 'splsda'
data(srbct)
X <- srbct$gene
Y <- srbct$class  
  
error <- valid(X, Y, ncomp = 10, method = "splsda", M = 3, 
               keepX = rep(10, 10), pred.method = "all")

plot(error, type = "l")
}

\keyword{regression}
\keyword{multivariate}
\keyword{hplot}
