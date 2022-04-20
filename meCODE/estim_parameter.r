################################################################################
#################################estim_parameter.r##############################
################################################################################
# R-Script to estimate h2, R and sigma_SNP with given u,v
#
#-------------------------------------------
#USAGE:
#       estim_param(b,u,v)
#
#
#--------------------------------------------
#INPUT:
#   u.............vector of breeding values
#   v.............vecor of SNP-effects
#
#-----------------------------------------------------------------
#OUTPUT:
#   r........estimated value for r
#   h2_3.....estimated h2 (HYB)
#   h2_1.....estimated h2 (PHEN)
#   h2_2.....estimated h2 (PED)
#   sigHYB.....estimated sigHYB
#   sigSNP...estimated sigSNP
#
#
# Himmelbauer, 15.12.2021
##################################################################################
library(data.table)
library(ggplot2)
library(Matrix)

estim_param=function(x){
    b<-         x[1:nf]
    u<-         x[(nf+1):(nf+na)]
    v<-         x[(nf+na+1):(nf+na+nSNP)]
    r<-         x[(nf+na+nSNP)+1]
    h<-         x[(nf+na+nSNP)+2]
    sigHYB<-    x[(nf+na+nSNP)+3]
    sigSNP<-    x[(nf+na+nSNP)+4]
    
    vc<-        m/(m-n)           #corraction for variance components
    
    #estim R
    gen<-       ZS(v)             #genomic EBV
    pa<-        as.vector(P(u))   #PA-EBV
    resid<-     (u-gen)           #Resid
    w<-         rep(0, length(u)) #weights all 0
    w[A$as]<-   1                 #except for genotyped animals
    Restim<-    lm(resid~pa, weight=w)[[1]][2]   #regression coefficient=r 
   
    #estim h2_3 (HYB)
    Rd<-        rep(1, na)
    Rd[A$as]<-  Restim
    resid2<-    (u-gen-pa*Rd)
    nup<-       as.numeric(classDef[,5]==2)+as.numeric(classDef[,6]==2)
    pi<-        ((nup[A$ac]+2)/4)  
    h23<-        mean((resid2/(pi*sigHYB^2))^2)
    sigHYBestim<- mean((resid2/(pi*h^2))^2)
    h23c<-        h23*vc
    sigHYBestimc<-sigHYBestim*vc

    #estim h2_1 (PHE)
    resid3<-    y-rep(b,length(y))-u[A$ap]
    h21<-       1-mean(resid3^2)
    h21c<-      h21*vc
    
    #estim h2_2  (PED)
    resid4<-    u-pa
    h22<-       mean((resid4/pi)^2)
    h22c<-      h22*vc 
    
    #estim h2_4  (PED+HYB)
    h24<-        (h22*na+h23*nas)/(na+nas)
    h24c<-       (h22c*na+h23c*nas)/(na+nas)

    #estim sig_SNP
    sigSNPestim<-  sd(v)
    sigSNPestimc<- sigSNPestim*vc
    
    epa<-          c(Restim, h21, h22, h23, h24, sigHYBestim, sigSNPestim, h21c, h22c, h23c, h24c, sigHYBestimc, sigSNPestimc)
    return(epa)
}