# Copyright (C) 2009 
# Kim-Anh Lê Cao, French National Institute for Agricultural Research and 
# ARC Centre of Excellence ins Bioinformatics, Institute for Molecular Bioscience, University of Queensland, Australia
# Leigh Coonan, Student, University of Queensland, Australia
# Fangzhou Yao, Student, University of Queensland, Australia
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


print.pca <- function(x, ...) {
    
 
    per.var = as.vector(x$sdev/sum(x$sdev))
    cum.var=as.vector(cumsum(per.var))
    x$sdev=as.vector(x$sdev)
    names(x$sdev) = paste("PC", 1:length(x$sdev), sep = "")
    names(per.var) = paste("PC", 1:length(per.var), sep = "")
    names(cum.var) = paste("PC", 1:length(cum.var), sep = "")
    

    cat("Eigenvalues for the first ", x$ncomp, "principal components:", "\n")
    print(x$sdev[1:x$ncomp])
    cat("\n")
    
    
    if(!x$NA.X) {
        cat("Proportion of explained variance for the first ", x$ncomp, "principal components:", "\n")
        print(per.var[1:x$ncomp])
        cat("\n")

        cat("Cumulative proportion explained variance for the first ", x$ncomp, "principal components:", "\n")
        print(cum.var[1:x$ncomp])
        cat("\n")

        cat(" Other available components: \n", "-------------------- \n")
        cat(" reconstructed data: see object$rotation \n")
    }
}
