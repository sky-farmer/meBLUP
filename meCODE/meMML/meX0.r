################################################################################
#################################meX0.r#########################################
################################################################################
# R-Script to prepare initial value for x
#
#-------------------------------------------
#USAGE:
#       prepX0(option)
#
#
#--------------------------------------------
#INPUT:
#   option..........possible options
#         ...rand: random vector with correct mean and variance
#         ...tbv:  using tbvs (missing values are random)  
#         ...fix_h2_param_opt:  second step using solSNPs estimated using TBVs with optimal starting values for EBVs
#         ...fix_h2_param_rand: second step using solSNPs estimated using TBVs with random starting values for EBVs
#         ...EBV_1: random starting values for EBVs, possible for genomic and conventional evaluations
#         ...EBV_2: conventional EBVs as starting values for EBVs
#         ...EBV_22: conventional EBVs as starting values for EBVs and estimated SNP-effects as starting values for v
#-----------------------------------------------------------------
#OUTPUT:
#      x0....initial vector for x in optim()
#
#
# Himmelbauer, 15.12.2021
##################################################################################
library(data.table)
library(ggplot2)
library(Matrix)

prepX0=function(option){

    #-----------------------------------------------
    # initial values for x
    #-----------------------------------------------
    
    #------------------------------------------------------------------------------
    #Read Data
    tbv<-    fread(paste0(PATH, '.bv'))
    tbv<-    tbv[tbv$id%in%A$cal]
    tbv<-    tbv$true
    miss<-   na-length(tbv)
    tbv<-    c(rnorm(miss, mean(tbv),sd(tbv)), tbv)
    dat<-    fread(paste0(PATH, '.dat'))
    dat<-    dat$V3
    fix<-    abs(mean(dat))
    snp<-    rnorm(nSNP, mean(0), SIG_SNP)
    
    #------------------------------------------------------------------------------
    #Different options
    
    if (option=='tbv') {
        h2<-     H2  
        r<-      Rho
        sigHYB<- SIG_HYB
        sigSNP<- SIG_SNP
        x0<<-    c(fix, tbv, snp, r, sqrt(h2), sigHYB, sigSNP)
    }
    
    if (option=='rand') {
        fix<-    abs(mean(dat))
        ebv<-    rnorm(na,   mean(tbv), sd(tbv))
        snp<-    rnorm(nSNP, mean(snp), SIG_SNP)
        r<-      Rho
        h2<-     0.5
        sigHYB<- SIG_HYB
        sigSNP<- SIG_SNP
        x0<<-    c(fix, ebv, snp, r, sqrt(h2), sigHYB, sigSNP)
    }
    
    if (option=='EBV_1') {
        fix<-    abs(mean(dat))
        ebv<-    rnorm(na,   mean(dat), sd(dat))
        snp<-    rnorm(nSNP, mean(0), 0.01)
        r<-      0.000000001
        h2<-     H2  
        sigHYB<- SIG_HYB
        sigSNP<- SIG_SNP
        if (nSNP>0) {
        x0<<-    c(fix, ebv, snp, r, sqrt(h2), sigHYB, sigSNP)
        } else {
        x0<<-    c(fix, ebv, r, sqrt(h2), sigHYB)
        }
    }
    
    if (option=='EBV_2') { 
        fix<-    abs(fix_conv)
        ebv<-    EBVcal_conv
        snp<-    rnorm(nSNP, mean(0), 0.01)
        r<-      Rho
        h2<-     H2  
        sigHYB<- SIG_HYB
        sigSNP<- SIG_SNP
        if (nSNP>0) {
        x0<<-    c(fix, ebv, snp, r, sqrt(h2), sigHYB, sigSNP)
        } else {
        x0<<-    c(fix, ebv, r, sqrt(h2), sigHYB)
        }
    }
    
        if (option=='EBV_3') {
        fix<-    abs(fix_conv)
        ebv<-    EBVcal_conv
        snp<-    vITB
        r<-      r_estim
        h2<-     H2  
        sigHYB<- SIG_HYB
        sigSNP<- SIG_SNP
        if (nSNP>0) {
        x0<<-    c(fix, ebv, snp, r, sqrt(h2), sigHYB, sigSNP)
        } else {
        x0<<-    c(fix, ebv, r, sqrt(h2), sigHYB)
        }
    }
    
    if (option=='fix_h2_param_opt') { #V3 optimaler start
        fix<-    abs(mean(tbv))
        snp<-    fread(paste0(PATH, suffix,'.solSNP'), header=T); setkey(snp,'SNP')
        snp<-    snp$Mean
        r<-      Rho
        h2<-     H2  
        sigHYB<- SIG_HYB    
        sigSNP<- SIG_SNP
        x0<<-    c(fix, tbv, snp, r, sqrt(h2), sigHYB, sigSNP)
    }
    
    if (option=='fix_h2_param_rand') { #V3 random start
        fix<-    abs(mean(tbv))
        ebv<-    rnorm(na,   mean(tbv), sd(tbv))
        snp<-    fread(paste0(PATH, suffix,'.solSNP'), header=T); setkey(snp,'SNP')
        snp<-    snp$Mean
        snp<-    rnorm(nSNP, mean(snp), sd(snp))
        r<-      Rho 
        h2<-     H2  
        sigHYB<- SIG_HYB    
        sigSNP<- SIG_SNP
        x0<<-    c(fix, ebv, snp, r, sqrt(h2), sigHYB, sigSNP)
    }
    
    return(x0)

}       