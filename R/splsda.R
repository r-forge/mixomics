# Copyright (C) 2009 
# S�bastien D�jean, Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
# Ignacio Gonz�lez, Genopole Toulouse Midi-Pyrenees, France
# Kim-Anh L� Cao, French National Institute for Agricultural Research and 
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


splsda <-
function(X, 
         Y,		
         ncomp = 2, 
	 keepX = rep(ncol(X), ncomp),
         max.iter = 500,		 
         tol = 1e-06,
         ...)
{
    X = as.matrix(X)
	
    #-- validation des arguments --#
    if (length(dim(X)) != 2 || !is.numeric(X)) 
        stop("'X' must be a numeric matrix.")
     
    if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0)
        stop("invalid number of variates, 'ncomp'.")
		
    # / Testing the input Y
    if (is.null(dim(Y))) {
        Y = as.factor(Y)	
        ind.mat = unmap(as.numeric(Y))					
    }
    else {
        stop("'Y' should be a factor or a class vector.")						
    }		
    # \ Testing input Y
	
    n = nrow(X)
     
    if ((n != nrow(ind.mat))) 
        stop("unequal number of rows in 'X' and 'Y'.")

    result = spls(X, ind.mat, ncomp = ncomp, mode = "regression", keepX = keepX, 
                  max.iter = max.iter, tol = tol, ...)
       
    cl = match.call()
    cl[[1]] = as.name('splsda')
    result$call = cl
	 
    result$ind.mat = ind.mat
    result$names$Y = levels(Y)
    class(result) = "splsda"
    return(invisible(result))	
}

