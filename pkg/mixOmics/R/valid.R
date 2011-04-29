# Copyright (C) 2009 
# Sébastien Déjean, Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
# Ignacio González, Genopole Toulouse Midi-Pyrenees, France
# Kim-Anh Lê Cao, French National Institute for Agricultural Research and 
# Queensland Facility for Advanced Bioinformatics, University of Queensland, Australia
# Pierre Monget, Ecole d'Ingenieur du CESI, Angouleme, France
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.


#------------------------------------------------------#
#-- Includes valid for PLS, sPLS, PLS-DA and sPLS-DA --#
#------------------------------------------------------#

valid <-
function(X, 
         Y, 
         ncomp = min(6, ncol(X)), 
         method = c("pls", "spls", "plsda", "splsda"),
         pred.method = c("all", "max.dist", "class.dist", "centroids.dist", "mahalanobis.dist"),
         mode = c("regression", "invariant", "classic"),
         keepX = NULL, keepY = NULL, 
         validation = c("loo", "Mfold"),
         M = if(validation == "Mfold") 10 else nrow(X),
         max.iter = 500, 
         tol = 1e-06)
{

    method = match.arg(method)
		
    #--------------------- PLS and sPLS ---------------------#
    if (any(c("pls", "spls") == method)) {
    	
        #-- validation des arguments --#
        #-- do warning for mode + other warnings --#
        if (length(dim(X)) != 2) 
            stop("'X' must be a numeric matrix.")
			
        mode = match.arg(mode)			
        if ((method == 'spls') & (mode == 'invariant')) 
            stop("No 'invariant' mode with sPLS.")  
			
        validation = match.arg(validation)
         
        X = as.matrix(X)
        Y = as.matrix(Y)
         	
        n = nrow(X)
        p = ncol(X)
        q = ncol(Y)
         
        if (!is.numeric(X) || !is.numeric(Y)) 
            stop("'X' and/or 'Y' must be a numeric matrix.")
        	
        if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0 || ncomp > p)
            stop("Invalid number of components, 'ncomp'.")
         
        if (any(is.na(X)) || any(is.na(Y))) 
            stop("Missing data in 'X' and/or 'Y'. Use 'nipals' for dealing with NAs.")
			
        if(is.null(keepX) && method == "spls") keepX = c(rep(ncol(X), ncomp))
        if(is.null(keepY) && method == "spls") keepY = c(rep(ncol(Y), ncomp))
          
        #-- criteria --#
        press.mat = Ypred = array(0, c(n, q, ncomp))
        MSEP = R2 = matrix(0, nrow = q, ncol = ncomp) 
         
        #-- M fold or loo cross validation --#
        ##- define the folds
        if (validation == "Mfold") { 
            fold = split(sample(1:n), rep(1:M, length = n)) 
        } 
        else { 
            fold = split(1:n, rep(1:n, length = n)) 
            M = n
        }
         
        for (i in 1:M) {
            omit = fold[[i]]
            X.train = X[-omit, ]
            Y.train = Y[-omit, ]
            X.test = matrix(X[omit, ], nrow = length(omit))
		    Y.test = matrix(Y[omit, ], nrow = length(omit))
	        	
            X.train = scale(X.train, center = TRUE, scale = FALSE)
            xmns = attr(X.train, "scaled:center")
             
            Y.train = scale(Y.train, center = TRUE, scale = FALSE)
            ymns = attr(Y.train, "scaled:center")
             
            X.test = scale(X.test, center = xmns, scale = FALSE)
	        	
            #-- pls or spls --#
            if (method == "pls") {
                object = pls(X = X.train , Y = Y.train, ncomp = ncomp, 
                             mode = mode, max.iter = max.iter, tol = tol)
            } 
            else {
                object = spls(X = X.train , Y = Y.train, ncomp = ncomp, 
                              mode = mode, max.iter = max.iter, tol = tol, 
                              keepX = keepX, keepY = keepY)
            }
		     
            Y.hat = predict(object, X.test)$predict
             
            for (h in 1:ncomp) {
			    Y.mat = matrix(Y.hat[, , h], nrow = dim(Y.hat)[1], ncol= dim(Y.hat)[2])		
                Y.hat[, , h] = sweep(Y.mat, 2, ymns, FUN = "+")		
                press.mat[omit, , h] = (Y.test - Y.hat[, , h])^2
                Ypred[omit, , h] = Y.hat[, , h]
            }
           
        }  #end i	
         	
        #-- compute MSEP and/or RMSEP --#
        for (h in 1:ncomp) { 
            MSEP[, h] = apply(as.matrix(press.mat[, , h]), 2, mean, na.rm = TRUE)
            R2[, h] = diag(cor(Y, Ypred[, , h], use = "pairwise"))		
        }
        	
	    RMSEP = sqrt(MSEP)
        	
        colnames(MSEP) = colnames(RMSEP)= colnames(R2) = paste('ncomp', c(1:ncomp), sep = " ")
        rownames(MSEP) = rownames(RMSEP) = rownames(R2) = colnames(Y)        		
         
        #-- valeurs sortantes --#
        res = list(msep = MSEP, rmsep = RMSEP, r2 = R2)
     
    }

    #-------------------------- for plsda and splsda ------------------------#
    if (any(c("plsda", "splsda") == method)) {
    	
        #-- validation des arguments --#
        #-- do warning for mode + other warnings --#
        X = as.matrix(X)
		 
        if (length(dim(X)) != 2 || !is.numeric(X)) 
            stop("'X' must be a numeric matrix.")
	    	
        if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0)
            stop("invalid number of variates, 'ncomp'.")
			
        if(is.null(keepX) && method == "splsda") keepX = c(rep(ncol(X), ncomp))
	     	
        if (is.null(dim(Y))) {
            Y = as.factor(Y)					
        }
        else {
            stop("'Y' should be a factor or a class vector.")						
        }		
         
        pred.method = match.arg(pred.method, several.ok = TRUE)		
        if (any(pred.method == "all")) nmthd = 4 
        else nmthd = length(pred.method)
		
        cl = split(1:length(Y), Y)
        M = min(c(M, sapply(cl, length)))
        cl.fold = fold = omit = list()
     	
        error.mat = array(0, dim = c(ncomp, nmthd, M))
        error.fun = function(x, y) {
            error.vec = sweep(x, 1, y, FUN = "-")
            error.vec = (error.vec != 0)		
            error.vec = apply(error.vec, 2, sum)/length(y)
            return(error.vec)
        }
         
        for (k in 1:length(cl)) {
            fold[[k]] = split(sample(cl[[k]]), rep(1:M, length = length(cl[[k]]))) 
        }
         	
        for (k in 1:length(cl)) {
            id = sample(1:M)
            fold[[k]] = fold[[k]][id]  
        }
          
        for (k in 1:M) {
            for (i in 1:length(cl)) {
                cl.fold[[i]] = fold[[i]][[k]]
            }
            omit[[k]] = unlist(cl.fold)
        }
         
        for (i in 1:M) {	
            X.train = X[-omit[[i]], ]
            Y.train = Y[-omit[[i]]]
            X.test = matrix(X[omit[[i]], ], nrow = length(omit[[i]]))
     	    	
            X.train = scale(X.train, center = TRUE, scale = FALSE)
            xmns = attr(X.train, "scaled:center")
             
            X.test = scale(X.test, center = xmns, scale = FALSE)
	        	
            #-- plsda or splsda --#
            if (method == "plsda") {
                object = plsda(X = X.train, Y = Y.train, ncomp = ncomp, 
                               max.iter = max.iter, tol = tol)
            } 
            else {
                object = splsda(X = X.train, Y = Y.train, ncomp = ncomp, 
                                max.iter = max.iter, tol = tol, 
                                keepX = keepX)
            }
      	    	
            Y.hat = predict(object, X.test, method = pred.method)$class		
            error.mat[, , i] = sapply(Y.hat, error.fun, y = as.numeric(Y[omit[[i]]]))
        }  #end i	
    	 
        #-- compute the error --#
        res = apply(error.mat, 1:2, mean)
        rownames(res) = paste('ncomp', 1:ncomp, sep = " ")
        colnames(res) = names(Y.hat)		
    } 
     
    method = paste(method, "mthd", sep = ".")	
    class(res) = c("valid", method)
    return(invisible(res))	
}
