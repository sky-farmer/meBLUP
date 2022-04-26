################################################################################
#################################mePREPARE.r####################################
################################################################################
# R-Script/function to set parameters for MMML
#
#-------------------------------------------
#USAGE:
# mePREP(CONFIG)
#   CONFIG....configuration parameters
#
#--------------------------------------------
#
# Himmelbauer, 14.11.2021
##################################################################################
library(data.table)
library(ggplot2)
library(Matrix)
rnd <- function(x) trunc(x+sign(x)*0.5)


mePREP=function(config){
    
    nSNP<<- config$nSNP
    prt<-   config$prt
    na<<-   config$na<-    length(A$ac)  #number of animals
    nap<<-  config$nap<-   nrow(A$X)     #number of animals with phenotypes
    nf<<-   config$nf<-    ncol(A$X)     #number of fixed effects 
    nas<<-  config$nas<-   nrow(A$S)     #number of animals with genotypes
    ns<<-   config$ns<-    ncol(A$S)     #number of SNPs
    m<<-    config$m<-     (nap+na+nSNP)
    n<<-    config$n<-     (nf+na+nSNP) 
    q<-     config$q<-     (m-n)/(2*m)
    
    
    # Adjust matrix to the number of SNPs
    ##########################################
    if (nSNP>0) {
    if (nSNP>ns || nSNP<1) {
        print(print0("Number of SNPs in Datafile is ", ns,". Number of SNPs given in config is ", nSNP, "."))
        stop("Incorrect number of SNPs!")
    }
    sel<<-        rnd((0.5:(nSNP-0.5))*(ns/nSNP))
    S<<-          A$S[,sel]
    A$Scan<<-     A$Scan[,sel]
    } else {
    sel<<-        NULL
    S<<-          NULL
    A$Scan<<-     NULL
    }
    
    ##########################################
    #right hand side
    ##########################################
    Fy<<-         c(y, rep(0,na+nSNP))
    
    return(config)
}















