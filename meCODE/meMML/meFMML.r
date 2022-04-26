################################################################################
#################################meFMML.r#######################################
################################################################################
# R-Script/function to define f_mml
#
#-------------------------------------------
#USAGE:
# fmml(x)
#   x....point at which the function should be calculated
#
#--------------------------------------------
#
# Himmelbauer, 17.05.2021 (nach Neumaier 2021, meBLUP.m)
##################################################################################

library(data.table)
library(Matrix)
rnd <- function(x) trunc(x+sign(x)*0.5)
 
#MML cost function
fmml<-function(x, iter='YES'){
    verbose=1
    
    if (!iter=='NO') {
    #Save current solution
    x<-      x*scl
    param[fmmlcount,names(param):=as.list(c(REP, x))]
    }
    
    rx<-     r(x)
    sx<-     Cix(x)%*%rx
    rxsx<-   crossprod(rx,sx)
    sdlcCx<- sdlcC(x)
    
    f<-      as.vector(0.5*(rxsx)+q*2*sdlcCx)
    if (is.nan(f)) {
        f<-  Inf
    }
    f<-      min(10^15,f)
    
    if (verbose) {
        fwrite(data.table(it=fmmlcount, FMML=f), 'FMML.out', col.names=F, row.names=F,quote=F,sep=' ',dec='.', append=TRUE)
        fmmlcount<<-fmmlcount+1
    }
    
    return(f)
}
    



