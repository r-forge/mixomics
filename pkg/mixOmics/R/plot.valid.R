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


plot.valid <- 
function (x, 
          criterion = c("MSEP", "RMSEP", "R2"),
          pred.method = "all",
          xlab = "number of components", 
          ylab = NULL,
          cTicks = NULL,
          layout = NULL,		  
          ...)
{

    #-- plot for pls and spls ----------------------------------------#
    if (any(class(x) == "pls.mthd") | any(class(x) == "spls.mthd")) {
         
        if (!any(criterion %in% c("MSEP", "RMSEP", "R2"))) 
            stop("Choose a validation criterion: MSEP, RMSEP or R2.")
        y = switch(criterion, MSEP = x$msep, RMSEP = x$rmsep, R2 = x$r2)
        	
        if (is.null(ylab))
            ylab = switch(criterion, MSEP = "MSEP", RMSEP = "RMSEP", R2 = expression(R^~2))
         	
        nResp = nrow(y)  # Number of response variables
        nComp = ncol(y)  # Number of components
         
        if (nResp > 1) {
            if (is.null(layout)) {
                nRows = min(c(3, nResp))
                nCols = min(c(3, ceiling(nResp / nRows)))
                layout = c(nRows, nCols)
            }
            else {
                if (length(layout) != 2 || !is.numeric(layout) || any(is.na(layout)))
                    stop("'layout' must be a numeric vector of length 2.")
                nRows = layout[1]
                nCols = layout[2]		
           }
    		
            if (nRows * nCols < nResp) devAskNewPage(TRUE) 
            ynames = rownames(y)
        }
        else {
            ynames = "Y"		
        }
         
        val = comps = vector("numeric")
        varName = vector("character")
         	
        for (i in 1:nResp) {
            val = c(val, y[i, ])
            comps = c(comps, 1:nComp)
            varName = c(varName, rep(ynames[i], nComp))
        }
         
        df = data.frame(val = val, comps = comps, varName = varName)
        if (is.null(cTicks)) cTicks = 1:ncol(x[[1]])
        yList = list(relation = "free")
			
	} # end plot for pls and spls	
	
    #-- plot for plsda and splsda ----------------------------------------#
    if (any(class(x) == "plsda.mthd") | any(class(x) == "splsda.mthd")) {
     
        if (any(pred.method == "all")) pred.method = colnames(x) 
		 
        if (!any(pred.method %in% colnames(x))) 
            stop("Choose the prediction methods.")
			
        x = matrix(x[, pred.method], ncol = length(pred.method))
        	
        if (is.null(ylab))
            ylab = "error rate"
         	
        nResp = ncol(x)  # Number of prediction methods
        nComp = nrow(x)  # Number of components
         
        if (nResp > 1) {
            if (is.null(layout)) {
                nRows = min(c(2, nResp))
                nCols = min(c(2, ceiling(nResp / nRows)))
                layout = c(nRows, nCols)
            }
            else {
                if (length(layout) != 2 || !is.numeric(layout) || any(is.na(layout)))
                    stop("'layout' must be a numeric vector of length 2.")
                nRows = layout[1]
                nCols = layout[2]		
           }
 		 
            if (nRows * nCols < nResp) devAskNewPage(TRUE) 
        }
         
        ynames = pred.method         
        val = comps = vector("numeric")
        varName = vector("character")
          	
        for (i in 1:nResp) {
            val = c(val, x[, i])
            comps = c(comps, 1:nComp)
            varName = c(varName, rep(ynames[i], nComp))
        }
         
        df = data.frame(val = val, comps = comps, varName = varName)
        if (is.null(cTicks)) cTicks = 1:nComp
        yList = list()		 
			
    } # end plot for plsda and splsda	
	
    xyplot(val ~ comps | varName, data = df, xlab = xlab, ylab = ylab,	
	    scales = list(y = yList, x = list(at = cTicks)), 
        as.table = TRUE, layout = layout, ...)
}
