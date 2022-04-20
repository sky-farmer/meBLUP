################################################################################
#################################meSTART_EBV.r################################
################################################################################
# R-Script to solve meBLUP using MML
#
#-------------------------------------------
#USAGE:
#       R CMD BATCH '--args "PATH" "H2" "NSNP" '                meSTART_EBV.r
#  or   R CMD BATCH '--args "PATH" "H2" "NSNP" "RHO" "SIG_SNP"' meSTART_EBV.r
# 
#--------------------------------------------
#INPUT:
#   PATH..............full path to dataset and name of the dataset
#   H2................known true H2
#   NSNP..............number of SNPs used
#   ------------ and optional ------------------------
#   RHO...............starting value for Rho (default: 0.2)
#   SIG_SNP...........starting value for SIG_SNP (default: 0.01)
#
#-----------------------------------------------------------------
# Detailed descrition of input files 
#
#   DATASET.ped........pedigree file [ID SIRE DAM GENERATION] (missing parent=0)  
#   DATASET.dat........phenotypic data file [ID FIX PHENO]
#   DATASET.snp........genotype data file [ID GENOTYPES(=012 coded with blanks in beetween)]
#   DATASET.bv.........file with TBVs and EBVs from SingleStep [ID GENERATION TBV ssGBLUP-GEBVS]
#
#OUTPUT:
#   PATH.solani........file with IDs for training animals and their EBVs [ID EBV]
#   PATH.solani_test...file with IDs for test animals and their EBVs [ID EBV]
#   PATH.solsnp........file with all snp effects [SNPEFFECT]
#   PATH.solfix........file with all fixed effects [FIXED EFFECTS]
#   PATH.solpar........file with parameters [rho h2 sigSNP sigHYB]
#   
# Himmelbauer, 29.12.2021 | last update 20.04.2022
##################################################################################

PATH_meBLUP<<-       'meMML/'   # Path to scripts with functions for meBLUP (based on MML)

#Package check
## packages needed in the Script
packages = c("data.table", "Matrix", "reshape2", "ggplot2", "optimx", "runner")

## Load or install&load all packeges
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

#Load and define R functions
source(paste0(PATH_meBLUP,"meRead.r"))
source(paste0(PATH_meBLUP,"mePREPARE.r"))
source(paste0(PATH_meBLUP,"meFMML.r"))
source(paste0(PATH_meBLUP,"meGRAD.r"))
source(paste0(PATH_meBLUP,"grad_check.r"))
source(paste0(PATH_meBLUP,"estim_parameter.r"))
source(paste0(PATH_meBLUP,"plot_out.r"))
source(paste0(PATH_meBLUP,"meX0.r"))
source(paste0(PATH_meBLUP,"plot_corCands.r"))

#####################################################################################
# INPUT SECTION
args=        commandArgs(trailingOnly=TRUE)
if (length(args)<3 & length(args)>0) { stop('At least 3 Inputs PATH, H2 and NSNP should be given! (Optional Rho and SIG_SNP)')}
if (length(args)==3) { print('No Inputs for Rho and SIG_SNP given. Default will be used!') }
if (length(args)>3 & length(args)<5) { stop('If the optional inputs are used, all 3 options Rho and SIG_SNP should be given!') }

PATH=        as.character(args[1])
H2=          as.numeric(args[2])
NSNP=        as.numeric(args[3])
if (length(args)>2) {
    Rho=     as.numeric(args[4])
    SIG_SNP= as.numeric(args[5])   
} else {
    Rho=     0.2
    SIG_SNP= 0.01
}

#or set manually
#PATH<-              '/data/.../.../test'
#H2<-                0.3
#NSNP<-              'all' #or number (30000)  
#Rho<-               0.2
#SIG_HYB<-           1
#SIG_SNP<-           0.01


# Default parameters
SIG_HYB<<-           1          # Value for sigHYB is set to 1 in this appraoch
PRT<<-               1          # Print level (0 or 1; 1...detailed outputs are written)
MAXIT1<<-            6000       # maximum number of iterations for calculation of conv. EBVs
MAXIT2<<-            1000       # maximum number of iterations for the first optimization with genomic information
MAXIT3<<-            100        # maximum number of iterations for the second optimization with genomic information
MISS<<-              0          # Missing value in phenotype data
SUFFIX<<-            ''         # Suffix used in output files
METHOD<<-            'L-BFGS-B' # Soving method used in optim(): 'L-BFGS-B'  'CG' 'spg' ...
option<<-            'rand'     # option..........possible options
                                #         ...rand: random vector with correct mean and variance
                                #         ...tbv:  using tbvs (missing values are random)  
                                #         ...fix_h2_param_opt:  second step using solSNPs estimated using TBVs with optimal starting values for EBVs
                                #         ...fix_h2_param_rand: second step using solSNPs estimated using TBVs with random starting values for EBVs
                                #         ...EBV_1: random starting values for EBVs, possible for genomic and conventional evaluations
                                #         ...EBV_2: conventional EBVs as starting values for EBVs
                                #         ...EBV_22: conventional EBVs as starting values for EBVs and estimated SNP-effects as starting values for v
#####################################################################################


print(Sys.time())

######################################################################################################################
# Step1: Calculate conv. EBVs to get starting values
######################################################################################################################
NSNP1<<-             0 #conv EBVs without genomic information
suffix<-             '_EBV'  
option<-             'EBV_1' 
METHOD<-             'L-BFGS-B' 

#####################################################################################
config<-             list()
config$prt<-         PRT       # print level

#Define number of SNPs to be used (default is all given SNPs are used)
if (!NSNP1=='all') {
    config$nSNP<-    NSNP1
} else {
    config$nSNP<-    ncol(A$S)
}

#-----------------------------------------------
# meRead-Read data
#-----------------------------------------------
    tmp<-            meRead(PATH, PRT)
    y<<-             tmp[[1]] 
    A<<-             tmp[[2]]
    classDef<<-      tmp[[3]]
    if (length(tmp)>3) {
        allBV<<-     tmp[[4]]
    }

#-----------------------------------------------
# mePrepare-Prepare optimization
#-----------------------------------------------
    #Set parameter
    config<-         mePREP(config)

    #Define operators
    source(paste0(PATH_meBLUP,"meOPERATORS.r"))
    
#-----------------------------------------------
# Prepare initial vecotr for x
#-----------------------------------------------
    x0<-              prepX0(option)
    x0[nf+na+nSNP+3]<-SIG_HYB
#-----------------------------------------------
# Scale x0
#-----------------------------------------------
    scl<<-            meSCALE(x0)
    x1<-              x0/scl
    
#-----------------------------------------------
# meBLUP-Gradcheck
#-----------------------------------------------
    #fmmlcount<-1
    #gradcount<-1
    #gradcheck(fmml, grad, x1, 1, c(1,(1+na/2), (1+na+nSNP/2), (1+na+nSNP+1), (1+na+nSNP+2), (1+na+nSNP+3), (1+na+nSNP+4)))

    #--------------
    #Gradient to compare (only for small models)
    #--------------
    #library(rootSolve)
    #g<-             gradient(fmml, x0)
    
#-----------------------------------------------
# Solve with optim
#-----------------------------------------------      
    # x1....Initial values for the parameters to be optimized over
    # fmml..A function to be minimized (or maximized), with first argument the vector of parameters over which minimization is to take place. It should return a scalar result.
    # grad..A function to return the gradient for the "BFGS", "CG" and "L-BFGS-B" methods. If it is NULL, a finite-difference approximation will be used.
    # lb....lower bound
    # ub....upper bound
    fmmlcount<-1
    gradcount<-1
    REP<-      1
    
    #lower and upper bound
    if (nSNP>0) {
    lb<-                 c(rep(-10, nf), rep(-10, na), rep(-10,nSNP), 0, 0.1 , 0.01 , 0.0001)/scl
    ub<-                 c(rep( 10, nf), rep( 10, na), rep( 10,nSNP), 0.99, 0.99, 1    , 15)/scl
    } else {
    lb<-                 c(rep(-10, nf), rep(-10, na), rep(-10,nSNP), 0, 0.1 , 0.01)/scl
    ub<-                 c(rep( 10, nf), rep( 10, na), rep( 10,nSNP), 0.99, 0.99, 1)/scl
    }
      
    #Save solution for all iterations
    param<<-             as.data.table(matrix(numeric(), (MAXIT1*5) ,(length(x0)+1)))
    if (nSNP>0) {
    names(param)<-       c('Rep',paste0('fix',rep(1:nf)),paste0('u',rep(1:na)),paste0('snp',rep(1:nSNP)), 'r', 'h2', 'sigHYB', 'sigSNP' ) 
    } else {
    names(param)<-       c('Rep',paste0('fix',rep(1:nf)),paste0('u',rep(1:na)), 'r', 'h2', 'sigHYB')
    }
    
    #Optimization
    sol <-               optim(x1, fmml, grad_EBV, method=METHOD, lower = lb, upper = ub, control = list(maxit=MAXIT1))
    
    print(paste0(        'Iterationen: ',nrow(param[!is.na(param$Rep)])))
    print(Sys.time())
    
    #Output of optim
    print(sol$value)
    print(sol$counts)
    print(sol$convergence)
    print(sol$message) 
    
    #Correlations for training animals
    par<-                as.matrix(na.omit(param))
    tbv<-                fread(paste0(PATH,'.bv'))
    tbvcal<-             tbv$true[!tbv$gen==max(tbv$gen)]
    correl<-             NULL
    for (i in (1:(nrow(par)))) {
        xi<-             par[i,-1]
        u<-              xi[(nf+1):(nf+na)]
        u<-              tail(u, length(tbvcal))
        co<-             cor(tbvcal, u)
        correl<-         c(correl, co)
    }
    correl_1<-           data.table(Iteration=1:length(correl), Correlation=correl, variable='train')
    
    #Calculate EBV and correlations for candidates for all iterations
    EBVcan<-             NULL
    for (i in (1:(nrow(par)))) {
        xi<-             par[i,-1]
        uCan<-           candsBV(xi)
        if (i==1) { 
            EBVcan<-     uCan
        } else {
            EBVcan<-     cbind(EBVcan, uCan[,c('EBV')])
        }
    }   
    tbvcan<-             tbv$true[tbv$gen==max(tbv$gen)]
    correl<-             NULL
    for (i in (2:(ncol(EBVcan)))) {
        c<-              cor(tbvcan, EBVcan[,..i])
        correl<-         c(correl, c)
    }
    correl_2<-           data.table(Iteration=1:length(correl), Correlation=correl, variable='test')
    
    #Plot correlations for conv. EBVs
    correl<-             rbind(correl_1, correl_2)
    p <-                 ggplot(correl,aes(Iteration, Correlation, linetype=variable))+
                            geom_line(size=0.8)+
                            scale_linetype_manual(values=c('solid','dotted'), labels=c('correl_train','correl_test'),name='')+
                            scale_size_manual(values=c(0.8,0.8),              labels=c('correl_train','correl_test'),name='')+
                            ylab('Correlation (Test and Train)')+
                            theme(axis.text=element_text(size=14),
                                  axis.title=element_text(size=16,face="bold"),
                                  legend.title = element_text(size=16), #change legend title font size
                                  legend.text = element_text(size=16),
                                  legend.position="bottom",
                                  legend.key.width = unit(1.9, "cm"),
                                  legend.key.height = unit(0.7, "cm"))
    png(paste0(PATH,'Correlation_konv',suffix,'.png'),width=12*300,height=7*300,res=300);plot(p);dev.off()
    
    
    #Save conv. EBVs
    xout<-               sol$par*scl
    fix_conv<<-          xout[(1):(nf)]
    EBVcal_conv<<-       xout[(nf+1):(nf+na)]
    EBVcan_conv<<-       candsBV(xout)
    tbv<-                fread(paste0(PATH,'.bv'))
    tbvcal<-             tbv$true[!tbv$gen==max(tbv$gen)]
    EBVconv<-            tail(EBVcal_conv, length(tbvcal))
    print(paste0('Correlation TBV and EBVs for training animals:',cor(tbvcal, EBVconv)))
    tbvcan<-             tbv$true[tbv$gen==max(tbv$gen)]
    print(paste0('Correlation TBV and EBVs for test animals:',cor(tbvcan, EBVcan_conv$EBV)))
    
    rm(param)
    
    
######################################################################################################################
# Step2: Using genomic information to get optimal SNP-Effects from fixed EBVs
######################################################################################################################
suffix<-             '_GEBV1' 
option<-             'EBV_2'
METHOD<-             'L-BFGS-B'
#####################################################################################
config<-             list()
config$prt<-         PRT       # print level

#-----------------------------------------------
# meRead-Read data
#-----------------------------------------------
    #Define number of SNPs to be used (default is all given SNPs are used)
    if (!NSNP=='all') {
        config$nSNP<- NSNP
    } else {
        config$nSNP<- ncol(A$S)
    }
    
    tmp<-             meRead(PATH, PRT)
    y<<-              tmp[[1]] 
    A<<-              tmp[[2]]
    classDef<<-       tmp[[3]]
    if (length(tmp)>3) {
        allBV<<-      tmp[[4]]
    }

#-----------------------------------------------
# mePrepare-Prepare optimization
#-----------------------------------------------
    #Set parameter
    config<-          mePREP(config)

    #Define operators
    source(paste0(PATH_meBLUP,"meOPERATORS.r"))
    
#-----------------------------------------------
# Prepare initial vecotr for x
#-----------------------------------------------
    x0<-              prepX0(option)
    x0[nf+na+nSNP+3]<-SIG_HYB
#-----------------------------------------------
# Scale x0
#-----------------------------------------------
    scl<<-            meSCALE(x0)
    x1<-              x0/scl
    
#-----------------------------------------------
# meBLUP-Gradcheck
#-----------------------------------------------
    #fmmlcount<-1
    #gradcount<-1
    #gradcheck(fmml, grad, x1, 1, c(1,(1+na/2), (1+na+nSNP/2), (1+na+nSNP+1), (1+na+nSNP+2), (1+na+nSNP+3), (1+na+nSNP+4)))

    #--------------
    #Gradient to compare (only for small models)
    #--------------
    #library(rootSolve)
    #g<-             gradient(fmml, x0)
    
#-----------------------------------------------
# Solve with optim
#-----------------------------------------------      
    # x1....Initial values for the parameters to be optimized over
    # fmml..A function to be minimized (or maximized), with first argument the vector of parameters over which minimization is to take place. It should return a scalar result.
    # grad..A function to return the gradient for the "BFGS", "CG" and "L-BFGS-B" methods. If it is NULL, a finite-difference approximation will be used.
    # lb....lower bound
    # ub....upper bound
    fmmlcount<-1
    gradcount<-1
    REP<-      1
    
     #lower and upper bound
    if (nSNP>0) {
    lb<-                 c(rep(-10, nf), rep(-10, na), rep(-10,nSNP), 0, 0.1 , 0.01 , 0.0001)/scl
    ub<-                 c(rep( 10, nf), rep( 10, na), rep( 10,nSNP), 1, 1,       1 , 15)/scl
    } else {
    lb<-                 c(rep(-10, nf), rep(-10, na), rep(-10,nSNP), 0, 0.1 , 0.01)/scl
    ub<-                 c(rep( 10, nf), rep( 10, na), rep( 10,nSNP), 0.99, 0.99, 1)/scl
    }
      
    #Save solution for all iterations
    param<<-             as.data.table(matrix(numeric(), (MAXIT2*5) ,(length(x0)+1)))
    if (nSNP>0) {
    names(param)<-       c('Rep',paste0('fix',rep(1:nf)),paste0('u',rep(1:na)),paste0('snp',rep(1:nSNP)), 'r', 'h2', 'sigHYB', 'sigSNP' ) 
    } else {
    names(param)<-       c('Rep',paste0('fix',rep(1:nf)),paste0('u',rep(1:na)), 'r', 'h2', 'sigHYB')
    }
    
    #Optimization
    sol <-               optim(x1, fmml, grad_GEBV1, method=METHOD, lower = lb, upper = ub, control = list(maxit=MAXIT2))
    
    print(paste0(        'Iterationen: ',nrow(param[!is.na(param$Rep)])))
    print(Sys.time())
    
    #Output of optim
    print(sol$value)
    print(sol$counts)
    print(sol$convergence)
    print(sol$message) 
    
    par<-                as.matrix(na.omit(param))
    #Calculate EBVs and correlations for candidates for all iterations
    EBVcan<-             NULL
    for (i in (1:(nrow(par)))) {
        xi<-             par[i,-1]
        uCan<-           candsBV(xi)
        if (i==1) { 
            EBVcan<-     uCan
        } else {
            EBVcan<-     cbind(EBVcan, uCan[,c('EBV')])
        }
    } 
    tbv<-                fread(paste0(PATH,'.bv'))
    print(paste0(        'Iterationen: ',(ncol(EBVcan)-1)))
    setkey(EBVcan,       'ID')
    setkey(tbv,          'id')
    tbv<-                tbv[J(EBVcan$ID)]
    correl<-             NULL
    for (i in (2:(ncol(EBVcan)))) {
        c<-              cor(tbv$true, EBVcan[,..i])
        correl<-         c(correl, c)
    }
    correl<-             data.table(Iteration=1:length(correl), Correlation=correl, variable='corr_test')    
        
    #Remove unused rows or param
    RES<-                NULL
    
    #Calculate EBV-PA for test animals
    for (i in (1:(nrow(par)))) {
        xi<-             par[i,-1]
        u<-              xi[(nf+1):(nf+na)]
        v<-              xi[(nf+na+1):(nf+na+nSNP)]
        r<-              xi[(nf+na+nSNP+1)]
        i<-              i+1
        ucan<-           as.vector(EBVcan[,..i]$EBV)
        Rdcan<-          data.table(ID=A$can, resid=rep(1, length(A$can)))
        Rdcan$resid[Rdcan$ID%in%A$can]<-r
        Rdcan<-          Rdcan$resid
        uPA<-            A$PAcan%*%u
        res<-            sd(ucan-uPA)
        RES<-            c(RES, res)
    }
    res<-                data.table(Iteration=1:length(res), value=RES, variable='sd(EBVtest-PA)')
    
    #Plot correlations and EBV-PA for test animals
    correl$logRES<-      log10(RES)
    ITB1<<-              correl$Iteration[correl$Correlation==max(correl$Correlation)][1]
    ITB2<<-              correl$Iteration[correl$logRES==min(correl$logRES)][1]
    ITB3<<-              round(nrow(par)/3)
    BEST<-               ITB1+1
    best<<-              EBVcan[,..BEST]
    correl2<-            NULL
    for (i in (2:(ncol(EBVcan)))) {
        c<-              cor(best, EBVcan[,..i])
        correl2<-        c(correl2, c)
    }
    correl2<-            data.table(Iteration=1:length(correl2), Correlation=correl2)
    ITB4<<-              correl2$Iteration[correl2$Correlation>=0.99][1]
    
    coeff<-              max(correl$logRES)*1.1
    
    #PLOT
    p<-                  ggplot(correl, aes(x=Iteration)) +
                            geom_line( aes(y=Correlation, linetype='Correlation'),    size=0.8) + 
                            geom_line( aes(y=logRES/coeff, linetype='log(sd(u-PA))'), size=0.8, ) +
                            geom_vline( aes(xintercept=ITB1,  color='optim_corr_test'),linetype='solid', size=0.5)+
                            geom_vline( aes(xintercept=ITB2,  color='optim_u-PA')     ,linetype='solid', size=0.5)+
                            geom_vline( aes(xintercept=ITB3,  color='Iter/3') ,linetype='solid', size=0.5)+
                            geom_vline( aes(xintercept=ITB4,  color='optim_corr_test0.99')     ,linetype='solid', size=0.5)+
                            scale_color_manual(name="", values=c("pink","red","green","blue"))+
                            scale_y_continuous(
                                name = "Correlation (Test)",
                                sec.axis = sec_axis(~.*coeff, name="log(sd(u-PA))")
                                ) + 
                            theme(axis.text=element_text(size=14),
                                  axis.title=element_text(size=16,face="bold"),
                                  legend.title = element_text(size=16),
                                  legend.text = element_text(size=16),
                                  legend.position="bottom",
                                  legend.key.width = unit(1.9, "cm"),
                                  legend.key.height = unit(0.7, "cm"),
                                  legend.box="vertical", legend.margin=margin())+
                            guides(color     = guide_legend(order = 0),
                                   linetype  = guide_legend(order = 1))+
                            scale_linetype_manual(values=c('solid','dashed'),name='')
    png(paste0(PATH,'RESIDCorPlot_Cands',suffix,'.png'),width=12*300,height=7*300,res=300);plot(p);dev.off()
    
    print(paste0('highest Correlation TBV-meBLUP at iteration:', correl[correl$Correlation==max(correl$Correlation),c('Iteration')]))
    print(paste0('highest Correlation TBV-meBLUP:',              correl[correl$Correlation==max(correl$Correlation),c('Correlation')]))
    
    #Plot sig_SNP for all iterations
    sigSNPhist<-            as.data.table(par[,c('sigSNP', 'Rep')])
    sigSNPhist$Iteration<-  1:nrow(sigSNPhist)
    sigSNPhist$sigSNP<-     log10(sigSNPhist$sigSNP^2)
    sigSNPhist<-            sigSNPhist[sigSNPhist$Iteration>4]
  
    sigSNPp <-              ggplot(sigSNPhist, aes(Iteration, sigSNP))+
                                geom_line(size=0.8) + 
                                geom_vline( aes(xintercept=ITB1,  color='optim_corr_test'),linetype='solid', size=0.5)+
                                geom_vline( aes(xintercept=ITB2,  color='optim_u-PA')     ,linetype='solid', size=0.5)+
                                geom_vline( aes(xintercept=ITB3,  color='Iter/3') ,linetype='solid', size=0.5)+
                                geom_vline( aes(xintercept=ITB4,  color='optim_corr_test0.99')     ,linetype='solid', size=0.5)+
                                scale_color_manual(name="", values=c("pink","red","green","blue"))+
                                ylab(expression(paste('log(',sigma^2,'_SNP)'))) +
                                theme(axis.text=element_text(size=14),  legend.title = element_text(size=16),
                                      legend.text = element_text(size=16),
                                      legend.position="bottom",
                                      legend.key.width = unit(1.9, "cm"),
                                      legend.key.height = unit(0.7, "cm"),
                                      axis.title=element_text(size=16,face="bold"),
                                      legend.box="vertical", legend.margin=margin())+
                                guides(color     = guide_legend(order = 0),
                                       linetype  = guide_legend(order = 1)) 
    png(paste0(PATH,'sigSNP_histPlot',suffix,'.png'),width=12*300,height=7*300,res=300);plot(sigSNPp);dev.off()
    
    #save solutions for iteration=niter/3
    ITb<-                   round(nrow(par)/3)
    uITB<<-                 par[ITb, (1+nf+1)   :(1+nf+na)]
    vITB<<-                 par[ITb, (1+nf+na+1):(1+nf+na+nSNP)]
    sSNPu<<-                sd(vITB)*1.5
    sSNPl<<-                sd(vITB)*0.5
    
    #Save a-posteriori estimated parameters
    tmp<-                   estim_param(par[ITb,-1])
    r_estim<<-              tmp[1]
    sigHYB_estim<-          tmp[6]
    sigSNP_estim<-          tmp[7]
    sigSNP_estimc<-         tmp[13]
    
######################################################################################################################
# Step3: fixed sigSNP, optimize u and v
######################################################################################################################
suffix<-              '_GEBV2'  
option<-              'EBV_3' 
METHOD<-              'L-BFGS-B'

#-----------------------------------------------
# mePrepare-Prepare optimization
#-----------------------------------------------
    #Set parameter
    config<-          mePREP(config)

    #Define operators
    source(paste0(PATH_meBLUP,"meOPERATORS.r"))
    
#-----------------------------------------------
# Prepare initial vecotr for x
#-----------------------------------------------
    x0<-              prepX0(option)
    x0[nf+na+nSNP+3]<-SIG_HYB
    
#-----------------------------------------------
# Scale x0
#-----------------------------------------------
    scl<<-            meSCALE(x0)
    x1<-              x0/scl
    
#-----------------------------------------------
# Solve with optim
#-----------------------------------------------      
    # x1....Initial values for the parameters to be optimized over
    # fmml..A function to be minimized (or maximized), with first argument the vector of parameters over which minimization is to take place. It should return a scalar result.
    # grad..A function to return the gradient for the "BFGS", "CG" and "L-BFGS-B" methods. If it is NULL, a finite-difference approximation will be used.
    # lb....lower bound
    # ub....upper bound
    fmmlcount<-1
    gradcount<-1
    REP<-      1
    
    #lower and upper bound
    if (nSNP>0) {
    lb<-               c(rep(-10, nf), rep(-10, na), rep(-10,nSNP), 0, 0.1 , 0.01 , sSNPl)/scl
    ub<-               c(rep( 10, nf), rep( 10, na), rep( 10,nSNP), 0.99, 0.99, 1    , sSNPu)/scl
    } else {
    lb<-               c(rep(-10, nf), rep(-10, na), rep(-10,nSNP), 0, 0.1 , 0.01)/scl
    ub<-               c(rep( 10, nf), rep( 10, na), rep( 10,nSNP), 0.99, 0.99, 1)/scl
    }
      
    #Save solution for all iterations
    param<<-           as.data.table(matrix(numeric(), (MAXIT1*5) ,(length(x0)+1)))
    if (nSNP>0) {
    names(param)<-     c('Rep',paste0('fix',rep(1:nf)),paste0('u',rep(1:na)),paste0('snp',rep(1:nSNP)), 'r', 'h2', 'sigHYB', 'sigSNP' ) 
    } else {
    names(param)<-     c('Rep',paste0('fix',rep(1:nf)),paste0('u',rep(1:na)), 'r', 'h2', 'sigHYB')
    }
    
    #Optimization
    sol <-             optim(x1, fmml, grad_GEBV2, method=METHOD, lower = lb, upper = ub, control = list(maxit=MAXIT1))
    
    print(paste0(      'Iterationen: ',nrow(param[!is.na(param$Rep)])))
    print(Sys.time())
    
    #Output of optim
    print(sol$value)
    print(sol$counts)
    print(sol$convergence)
    print(sol$message)     
    
#-----------------------------------------------
# Remove unused rows or param
#-----------------------------------------------
    param<-            na.omit(param)
    
#-----------------------------------------------
# Estim parameter for all iterations
#-----------------------------------------------
    estim<-            NULL
    par<-              as.matrix(param)
    for (i in (1:(nrow(param)))) {
        tmp<-          estim_param(par[i,-1])
        r_estim<-      tmp[1]
        h21<-          tmp[2]
        h22<-          tmp[3]
        h23<-          tmp[4]
        h24<-          tmp[5]
        sigHYB_estim<- tmp[6]
        sigSNP_estim<- tmp[7]
        h21c<-         tmp[8]
        h22c<-         tmp[9]
        h23c<-         tmp[10]
        h24c<-         tmp[11]
        sigHYB_estimc<-tmp[12]
        sigSNP_estimc<-tmp[13]
        e<-            data.table(estim_r=r_estim, h21=h21, h22=h22, h23=h23,h24=h24, h21c=h21c, h22c=h22c, h23c=h23c, h24c=h24c, estim_sigHYB=sigHYB_estim, estim_sigSNP=sigSNP_estim, estim_sigHYBc=sigHYB_estimc, estim_sigSNPc=sigSNP_estimc)
        estim<-        rbind(estim, e)
    }
    param<-            cbind(param, estim)
    
#-----------------------------------------------
# Calculate EBV for candidates for all iterations
#-----------------------------------------------
    EBVcan<-           NULL
    for (i in (1:(nrow(param)))) {
        xi<-           par[i,-1]
        uCan<-         candsBV(xi)
        if (i==1) { 
            EBVcan<-   uCan
        } else {
            EBVcan<-   cbind(EBVcan, uCan[,c('EBV')])
        }
    }   

#-----------------------------------------------
# Write output
#----------------------------------------------- 
    if (PRT==1) {    
        #Write output for all iterations
        fwrite(EBVcan,  paste0(PATH,'.solaniCanALL',suffix), col.names=TRUE, row.names=FALSE,quote=FALSE,sep=' ',dec='.') 
        #Write output for all iterations
        fwrite(as.data.table(param),   paste0(PATH,'.ebv_iter',suffix), col.names=TRUE, row.names=FALSE,quote=FALSE,sep=' ',dec='.')    
    }
    
    #Plot correlation for test animals
    plot_corCan(EBVcan, param, TRUE)
    
    #Plot output (EBVs fro training animals and parameters)
    plot_out(param, PATH, suffix)
    
    #Write output files solani, solsnp, solfix and solpar for best solution
    solpar<-           param[nrow(param), c('r', 'h2', 'sigSNP', 'sigHYB')]
    solfix<-           param[nrow(param), (1+1):nf]
    solani<-           data.table(ID=A$cal, ebv_m=t(param[nrow(param),(nf+1+1):(nf+na+1)]))
    solaniTest<-       data.table(ID=um_test$ID, ebv_m=EBVcan[,ncol(EBVcan)])
    solsnp<-           param[nrow(param),(nf+na+1+1):(nf+na+nSNP)]
    fwrite(solpar,     paste0(PATH,'.solpar',suffix),      col.names=TRUE,  row.names=FALSE,quote=FALSE,sep=' ',dec='.')
    fwrite(solfix,     paste0(PATH,'.solfix',suffix),      col.names=FALSE, row.names=FALSE,quote=FALSE,sep=' ',dec='.')
    fwrite(solani,     paste0(PATH,'.solani',suffix),      col.names=FALSE, row.names=FALSE,quote=FALSE,sep=' ',dec='.')
    fwrite(solaniTest, paste0(PATH,'.solani_test',suffix), col.names=FALSE, row.names=FALSE,quote=FALSE,sep=' ',dec='.')
    fwrite(solsnp,     paste0(PATH,'.solsnp',suffix),      col.names=FALSE, row.names=FALSE,quote=FALSE,sep=' ',dec='.')   

#-----------------------------------------------
# Compare GEBVs with TBVs and GEBVs from single step
#-----------------------------------------------    
    #GEBVS from meBLUP
    um<-               as.data.table(param)
    um<-               data.table(ID=A$cal, ebv_m=t(um[nrow(um),(nf+1+1):(nf+na+1)]))
    names(um)<-        c('ID', 'ebv_m')
    um_test<-          EBVcan; x<-ncol(um_test)
    um_test<-          data.table(ID=um_test$ID, ebv_m=um_test[,..x])
    names(um_test)<-   c('ID', 'ebv_m')
    
    #TBVs
    u<-                fread(paste0(PATH,'.bv'))
    u<-                data.table(ID=u$id, tbv=u$true, sstep=u$sstep)
    
    u_test<-           u[u$ID%in%um_test$ID]
    u<-                u[!u$ID%in%u_test$ID]
    u$meBLUP<-         um$ebv_m
    u_test$meBLUP<-    um_test$ebv_m

    #Pairwise correaltions (training)
    cor(u$tbv,         u$sstep)
    cor(u$tbv,         u$meBLUP)
    cor(u$sstep,       u$meBLUP)

    #Pairwise correaltions (test)
    cor(u_test$tbv,    u_test$sstep)
    cor(u_test$tbv,    u_test$meBLUP)
    cor(u_test$sstep,  u_test$meBLUP)
    
    
print(Sys.time())
warnings()
    
q('no')








