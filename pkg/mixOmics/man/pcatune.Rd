\name{pcatune}
\encoding{latin1}
\alias{pcatune}

\title{Tune the number of principal components in PCA}

\description{ 
'pcatune' can be used to quickly visualise the eigenvalues and the proportions of explained varience 
to choose the number of principal components to analyse in 'pca'. If the data is complete 'pcatune' 
uses Singular Value Decomposition and if there is data missing it uses the NIPALS algorithm. 
}

\usage{
pcatune(X, ncomp = NULL, center = TRUE,  scale. = FALSE,  max.iter = 500, tol = 1e-09)
}

\arguments{
  \item{X}{a numeric matrix (or data frame) which provides the 
    data for the principal components analysis. It can contain missing values.}
  \item{ncomp}{integer, the number of components to initially analyse in 'pcatune' to choose a final 'ncomp' for 'pca'. If \code{NULL}, 
     function sets ncomp = min(nrow(X),ncol(X))}	
  \item{center}{a logical value indicating whether the variables should be shifted to be zero centered. 
    Alternately, a vector of length equal the number of columns of \code{X} can be supplied. 
    The value is passed to \code{\link{scale}}.} 
  \item{scale.}{a logical value indicating whether the variables should be scaled to have 
    unit variance before the analysis takes place. The default is FALSE for consistency with \code{prcomp} 
    function, but in general scaling is advisable. Alternatively, a vector of length equal the number of 
    columns of \code{X} can be supplied. The value is passed to \code{\link{scale}}.} 
  \item{max.iter}{integer, the maximum number of iterations in the nipals algorithm.}
  \item{tol}{a positive real, the tolerance used in the nipals algorithm.}
}

\details{
The calculation is done either by a singular value decomposition of the (possibly centered and scaled) 
data matrix, if the data is complete or by using the nipals algorithm if there is data missing. Unlike 
\code{\link{princomp}}, the print method for these objects prints the results in a nice format and the 
\code{plot} method produces a bar plot of the percentage of variance explained by the principal 
components (PCs).

Note that \code{scale. = TRUE} cannot be used if there are zero or constant (for \code{center = TRUE}) variables. 
 
}    

\value{
\code{pcatune} returns a list with class "pcatune" containing the following components: 
  \item{eigenvalues}{The eigenvalues of the covariance/correlation matrix, though the calculation is 
actually  done with the singular values of the data matrix).} 
  \item{Proportion of explained variance}{the proportion of explained variance accounted for by each 
principal component is calculated using the eigen values}
  \item{Cumulative proportion of explained variance}{The cumulative proportion of explained variance 
accounted for by the sequential accumulation of principal components is calculated using the sum of the
 'proportion of explained variance'}
  \item{Estimated proportion of explained variance}{the estimated proportion of explained variance 
accounted for by each principal component is calculated using the first (min(ncol(X),nrow(X)) eigen 
values under the assumption that the first (min(ncol(X),nrow(X))principal components will account for
 100 \% of the explained variance}
  \item{Estimated cumulative proportion of explained variance}{The cumulative proportion of explained 
variance accounted for by the sequential accumulation of principal components is calculated using the 
sum of the 'proportion of explained variance'}
}

\author{Ignacio González and Leigh Coonan}

\seealso{\code{\link{nipals}}, \code{\link{biplot}}, 
\code{\link{plotIndiv}}, \code{\link{plotVar}},
\code{\link{plot3dIndiv}}, \code{\link{plot3dVar}}.}

\keyword{algebra}
