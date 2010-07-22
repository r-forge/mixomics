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


splsda <-
function(X, 
         Y,		
         ncomp = 2, 
		 keepX = c(rep(ncol(X), ncomp)),
         max.iter = 500,
         scaleY = TRUE,
         mode = "regression",		 
         tol = 1e-06)
{
	#-- validation des arguments --#
    if (length(dim(X)) != 2 || !is.numeric(X)) 
        stop("'X' must be a numeric matrix.")
     
    n = nrow(X)
     
    if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0)
        stop("invalid number of variates, 'ncomp'.")
		
	# / Testing of the inputs 		
	if(is.null(dim(Y))){			
			if(is.vector(Y)){
				if(is.integer(Y)){				
				Yprim = unmap(as.numeric(as.factor(Y)))						
				}else {stop("the vector 'Y' must be made of integers")
				}		
			}
			if(is.factor(Y)){
			Yprim = unmap(as.numeric(Y))			
			}	 		
	}else{
		if (is.matrix(Y)){ 
		stop("Y is a matrix, it should be a class vector !")
		}else{
		stop("You should enter either a vector or a factor !")
		}
	}	
	# \ Testing of the inputs

	result = spls(X, Yprim, ncomp = ncomp, mode = "regression", keepX = keepX, max.iter = max.iter, 
                 tol = tol, scaleY = scaleY)

	result$Yprim = Yprim
    class(result) = "pls"
    return(invisible(result))	
}

