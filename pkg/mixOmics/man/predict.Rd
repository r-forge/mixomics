\name{predict}
\encoding{latin1}
\alias{predict.splsda}
\alias{predict.plsda}

\title{Predict Method for PLSDA Regression and Sparse PLSDA}

\description{
Predicted values based on PLSDA regression or sparse PLSDA models. New responses and 
variates are predicted using a fitted model and a new matrix of observations.
}

\usage{

\method{predict}{splsda}(object, newdata, method = c("class.dist", "centroids.dist", "Sr.dist", "KA.dist"), ...)

\method{predict}{plsda}(object, newdata, method = c("class.dist", "centroids.dist", "Sr.dist", "KA.dist"), ...)
}

\arguments{
  \item{object}{object of class inheriting from \code{"splsda"} or \code{"plsda"}.}
  \item{newdata}{data matrix in which to look for for explanatory variables to be used for prediction.}
  \item{method}{method to aplly for \code{"splsda"} or \code{"plsda"}. It can be either "Ka.dist", "centroids.dist", "class.dist", or "Sr.dist".}
  \item{...}{other arguments to be passed \code{predict.default}.}
}

\details{
\code{predict} produces predicted values, obtained by evaluating the PLSDA/sparse PLSDA  
model returned by \code{plsda} or \code{splsda} in the frame \code{newdata}. 
Variates for \code{newdata} are also returned. 
}

\value{
\code{predict} produces a list with the following components: 
  \item{predict}{a three dimensional array of predicted response values. The dimensions 
	correspond to the observations, the response variables and the model dimension, respectively.}
  \item{variates}{matrix of predicted variates.}
  \item{B.hat}{matrix of regression coefficients (without the intercept).}
}

\references{
Tenenhaus, M. (1998). \emph{La régression PLS: théorie et pratique}. Paris: Editions Technic.  
}

\author{Sébastien Déjean, Ignacio González, Kim-Anh Lê Cao and Pierre Monget}

\seealso{\code{\link{plsda}}, \code{\link{splsda}}.}

\examples{
## First example
ncomp = 2 
max.iter = 500
tol = 1e-06

data(liver.toxicity)
X = as.matrix(liver.toxicity$gene)
Y = liver.toxicity$treatment[,4] 

i = 1
samp = sample(1:3, nrow(X), replace = TRUE)  
test = which(samp==i)
train = setdiff(1:nrow(X), test)

plsda.train = plsda(X[train,], Y[train], ncomp = ncomp, scaleY = TRUE,  mode = 'regression')
test.predict = predict(plsda.train, X[test,], method = "class.dist")
test.predict$cl

## Second example
ncomp = 2 
keepX = 20
max.iter = 500
tol = 1e-06

data(liver.toxicity)
X = as.matrix(liver.toxicity$gene)
Y = liver.toxicity$treatment[,4]  # time points

i = 1
samp = sample(1:3, nrow(X), replace = TRUE)  # here I perform 3-CV 
test = which(samp==i)
train = setdiff(1:nrow(X), test)

splsda.train = splsda(X[train,], Y[train], ncomp = ncomp, keepX = c(rep(keepX, ncomp)), scaleY = TRUE,  mode = 'regression')
test.predict = predict(splsda.train, X[test,], method = "class.dist")
test.predict$cl
}

\keyword{regression}
\keyword{multivariate}
