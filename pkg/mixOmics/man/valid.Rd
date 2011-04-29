\name{valid}
\encoding{latin1}
\alias{valid}

\title{Compute validation criterion for PLS, sPLS, PLS-DA and sPLS-DA}

\description{
Function to estimate measures of the prediction error for fitted PLS, sparse PLS,
PLS-DA and sparse PLS-DA models. M-fold and leave-one-out cross-validation are implemented.
}

\usage{
valid(X, Y, ncomp = min(6, ncol(X)), 
      method = c("pls", "spls", "plsda", "splsda"),
      pred.method = c("all", "max.dist", "class.dist", 
                      "centroids.dist", "mahalanobis.dist"),
      mode = c("regression", "invariant", "classic"),
      keepX = NULL, keepY = NULL, 
      validation = c("loo", "Mfold"),
      M = if(validation == "Mfold") 10 else nrow(X),
      max.iter = 500, 
      tol = 1e-06)
}

\arguments{
  \item{X}{numeric matrix of predictors.}
  \item{Y}{numeric vector or matrix of responses (for multi-response models).}
  \item{ncomp}{the number of components to include in the model. Default is \code{1:min(6, ncol(X))}.}
  \item{mode}{character string. What type of algorithm to use, matching one of \code{"classic"}, \code{"invariant"} 
    or \code{"regression"}.}
  \item{method}{Choose between \code{pls} or \code{spls}.}
  \item{keepX}{if \code{method="spls"} numeric vector of length \code{ncomp}, the number of variables
    weights to keep in \eqn{X}-loadings. By default all variables are kept in the model.}
  \item{keepY}{if \code{method="spls"} numeric vector of length \code{ncomp}, the number of variables
    weights to keep in \eqn{Y}-loadings. By default all variables are kept in the model.}
  \item{pred.method}{prediction method to be applied for \code{plsda} or \code{splsda}. Should be a subset 
    of \code{"max.dist"}, \code{"class.dist"}, \code{"centroids.dist"}, \code{"mahalanobis.dist"}. 
	Default is \code{"all"}. See \code{\link{predict}}.}
  \item{validation}{character.  What kind of (internal) validation to use.  See below.}
  \item{M}{the number of folds in the Mfold cross-validation.}
  \item{max.iter}{integer, the maximum number of iterations.}
  \item{tol}{a not negative real, the tolerance used in the iterative algorithm.}
}

\details{
For fitted PLS and sPLS regression models, \code{valid} estimates the 
mean squared error of prediction (MSEP), root mean squared error of prediction 
(RMSEP) and \eqn{R^2} to assess the predictive validity of the model using M-fold or 
leave-one-out cross-validation. Note that only the \code{classic}, \code{regression} and 
\code{invariant} modes can be applied.

If \code{validation = "Mfold"}, M-fold cross-validation is performed. 
How many folds to generate is selected by specifying the number of folds in \code{M}.
If \code{validation = "loo"}, leave-one-out cross-validation is performed.

For fitted PLS-DA and sPLS-DA models, \code{valid} estimates the classification error rate 
using cross-validation. How many folds to generate is selected such that there is at least 1 
sample for each class in the test set.
}

\value{
For PLS and sPLS models, \code{valid} produces a list with the following components: 
  \item{msep}{Mean Square Error Prediction for each Y variable.}
  \item{rmsep}{Root Mean Square Error Prediction for each Y variable.}
  \item{r2}{a matrix of \eqn{R^2} values of the \eqn{Y}-variables for models 
    with \eqn{1,...,}\code{ncomp} components.}
	
For PLS-DA and sPLS-DA models, \code{valid} produces a matrix. The dimensions 
	correspond to the components in the model and to the prediction method used, respectively.
}

\references{
Tenenhaus, M. (1998). \emph{La régression PLS: théorie et pratique}. Paris: Editions Technic.  

Lê Cao, K. A., Rossouw D., Robert-Granié, C. and Besse, P. (2008). A sparse PLS for variable 
selection when integrating Omics data. \emph{Statistical Applications in Genetics and Molecular 
Biology} \bold{7}, article 35.

Mevik, B.-H., Cederkvist, H. R. (2004). Mean Squared Error of Prediction (MSEP) Estimates for Principal Component 
Regression (PCR) and Partial Least Squares Regression (PLSR). \emph{Journal of Chemometrics} \bold{18}(9), 422-429.
}

\author{Sébastien Déjean, Ignacio González and Kim-Anh Lê Cao.}

\seealso{\code{\link{predict}}, \code{\link{nipals}}, \code{\link{plot.valid}}.}

\examples{
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
