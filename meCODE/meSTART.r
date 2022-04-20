################################################################################
#################################meSTART.r################################
################################################################################
# R-Script to solve meBLUP using MML
#
#-------------------------------------------
#USAGE:
#       R CMD BATCH '--args "PATH" "H2" "NSNP" '                          meSTART.r
#  or   R CMD BATCH '--args "PATH" "H2" "NSNP" "RHO" "SIG_HYB" "SIG_SNP"' meSTART.r
# 
#--------------------------------------------
#INPUT:
#   PATH..............full path to dataset and name of the dataset
#   H2................known true H2
#   NSNP..............number of SNPs used
#   --------------- and optional ------------------------
#   RHO...............starting value for Rho (default: 0.2)
#   SIG_HYB...........starting value for SIG_HYB (default: 1)
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
source(paste0(PATH_meBLUP,"plot_out_REP.r"))
source(paste0(PATH_meBLUP,"meX0.r"))
source(paste0(PATH_meBLUP,"plot_corCands.r"))
source(paste0(PATH_meBLUP,'plot_corCands_REP.r'))


#####################################################################################
# INPUT SECTION
args=        commandArgs(trailingOnly=TRUE)
if (length(args)<3 & length(args)>0) { stop('At least 3 Inputs PATH, H2 and NSNP should be given! (Optional Rho, SIG_HYB and SIG_SNP)')}
if (length(args)==3) { print('No Inputs for Rho, SIG_HYB and SIG_SNP given. Default will be used!') }
if (length(args)>3 & length(args)<6) { stop('If the optional inputs are used, all 3 options Rho, SIG_HYB and SIG_SNP should be given!') }

PATH=        as.character(args[1])
H2=          as.numeric(args[2])
NSNP=        as.numeric(args[3])
if (length(args)>2) {
    Rho=     as.numeric(args[4])
    SIG_HYB= as.numeric(args[5])
    SIG_SNP= as.numeric(args[6])   
} else {
    Rho=     0.2
    SIG_HYB= 1
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
PRT<<-               1          # Print level (0 or 1; 1...detailed outputs are written)
MAXIT1<<-            5000       # maximum number of iterations in the first optimization
MAXIT2<<-            500        # maximum number of iterations in repetitions
MISS<<-              0          # Missing value in phenotype data
SUFFIX<<-            ''         # Suffix used in output files
METHOD<<-            'L-BFGS-B' # Soving method used in optim(): 'L-BFGS-B'  'CG' 'spg' ...
REPMAX<<-            1          # Maximum number of repetitions
option<<-            'rand'     # option..........possible options
                                #         ...rand: random vector with correct mean and variance
                                #         ...tbv:  using tbvs (missing values are random)  
                                #         ...fix_h2_param_opt:  second step using solSNPs estimated using TBVs with optimal starting values for EBVs
                                #         ...fix_h2_param_rand: second step using solSNPs estimated using TBVs with random starting values for EBVs
                                #         ...EBV_1: random starting values for EBVs, possible for genomic and conventional evaluations
                                #         ...EBV_2: conventional EBVs as starting values for EBVs
                                #         ...EBV_22: conventional EBVs as starting values for EBVs and estimated SNP-effects as starting values for v
#####################################################################################

# Saving parameters
config<-             list()
config$maxit1<-      MAXIT1    # maximal number of iterations
config$maxit2<-      MAXIT2    # maximal number of iterations
config$prt<-         PRT       # print level

#-----------------------------------------------
# meRead-Read data
#-----------------------------------------------
    #Define number of SNPs to be used (default is all given SNPs are used)
    if (!NSNP=='all') {
        config$nSNP<-NSNP
    } else {
        config$nSNP<-ncol(A$S)
    }
    
    tmp<-            meRead(PATH, PRT)
    y<<-             tmp[[1]] 
    A<<-             tmp[[2]]
    classDef<<-      tmp[[3]]
    if (length(tmp)>3) {
        allBV<<-     tmp[[4]]
    }

#-----------------------------------------------
# mePrepare-Prepare optimization/Set parameters
#-----------------------------------------------
    config<-    mePREP(config)

#-----------------------------------------------  
#Define operators
#-----------------------------------------------
    source(paste0(PATH_meBLUP,"meOPERATORS.r"))
    
#-----------------------------------------------
# Prepare initial vecotr for x
#-----------------------------------------------
    x0<-        prepX0(option)
    
#-----------------------------------------------
# Scale x0
#-----------------------------------------------
    scl<<-      meSCALE(x0)
    x1<-        x0/scl
    
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
    lb<-                 c(rep(-10, nf), rep(-10, na), rep(-10,nSNP), 0.01, 0.1 , 0.01 , 0.0001)/scl
    ub<-                 c(rep( 10, nf), rep( 10, na), rep( 10,nSNP), 0.99, 0.99, 1    , 15)/scl
      
    #Save solution for all iterations
    param<<-         as.data.table(matrix(numeric(), ((MAXIT1+5*MAXIT2)*5) ,(length(x0)+1)))
    names(param)<-   c('Rep',paste0('fix',rep(1:nf)),paste0('u',rep(1:na)),paste0('snp',rep(1:nSNP)), 'r', 'h2', 'sigHYB', 'sigSNP' ) 
    
    #Optimization
    sol <-            optim(x1, fmml, grad, method=METHOD, lower = lb, upper = ub, control = list(maxit=config$maxit1))
    
    print(paste0(    'Iterationen: ',nrow(param[!is.na(param$Rep)])))
    print(Sys.time())
    
    #Output of optim
    print(sol$value)
    print(sol$counts)
    print(sol$convergence)
    print(sol$message)   

    #Solution
    xout<-            sol$par*scl
    er<-              xout[length(xout)-3]
    eh2<-             xout[length(xout)-2]^2
    esigHYB<-         xout[length(xout)-1]^2
    esigSNP<-         xout[length(xout)]^2
    print(paste0(     'h2=', eh2))
    print(paste0(     'r=', er))
    print(paste0(     'sigma_HYB=', esigHYB))
    print(paste0(     'sigma_SNP=', esigSNP))
    
    if (PRT==1) {
        #Write output for all iterations
        fwrite(as.data.table(param),   paste0(PATH,'.ebv_iter',suffix), col.names=TRUE, row.names=FALSE,quote=FALSE,sep=' ',dec='.')   
    }
#-----------------------------------------------
# Repeat the optimization with different starting values
#-----------------------------------------------
    while (REP<REPMAX) {
        REP<-       REP+1
    
        #--------------------------------
        #Search for the best solutions (highest corr_test)
        #--------------------------------
        par<-     na.omit(param)
        par<-     as.matrix(par)
        EBVcan<-  NULL
        for (i in (1:(nrow(par)))) {
            xi<-           par[i,-1]
            uCan<-         candsBV(xi)
            if (i==1) { 
                EBVcan<-   uCan
            } else {
                EBVcan<-   cbind(EBVcan, uCan[,c('EBV')])
            }
        }   
    
        tbv<-          fread(paste0(PATH,'.bv'))
        print(paste0(  'Iterationen: ',(ncol(EBVcan)-1)))
        setkey(EBVcan, 'ID')
        setkey(tbv,    'id')
        tbv<-          tbv[J(EBVcan$ID)]
        correl<-       NULL
    
        for (i in (2:(ncol(EBVcan)))) {
            c<-        cor(tbv$true, EBVcan[,..i])
            correl<-   c(correl, c)
        }
        correl<-    data.table(Iteration=1:length(correl), Correlation=correl, Rep=par[,1])
    
        print(paste0('highest Correlation TBV-meBLUP at iteration:', correl[correl$Correlation==max(correl$Correlation),c('Iteration')]))
        print(paste0('highest Correlation TBV-meBLUP:', correl[correl$Correlation==max(correl$Correlation),c('Correlation')]))
    
        ITbest<-    correl$Iteration[correl$Correlation==max(correl$Correlation)][1]
        xe<-        unlist(param[ITbest[length(ITbest)],-1])
    
        #--------------------------------
        #Schaetze Varianzkomponenten aus der besten Loesung
        #--------------------------------    
        tmp<-          estim_param(xe)
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
    
        print(paste0('r=',r_estim))
        print(paste0('h2_1(PHE)=',h21))
        print(paste0('h2_2(PED)=',h22))
        print(paste0('h2_3(HYB)=',h23))
        print(paste0('h2_4(PEd+HYB)=',h24))
        print(paste0('sigHYBestim=',sigHYB_estim))
        print(paste0('sigSNPestim=',sigSNP_estim))
        print(paste0('h2_1(PHE)*vc=',h21c))
        print(paste0('h2_2(PED)*vc=',h22c))
        print(paste0('h2_3(HYB)*vc=',h23c))
        print(paste0('h2_4(PED+HYB)*vc=',h24c))
        print(paste0('sigHYBestim*vc=',sigHYB_estimc))
        print(paste0('sigSNPestim*vc=',sigSNP_estimc))
      
        #Mean h2
        print(paste0('mean h2=',mean(c(h21, h22, h23, h21c, h22c, h23c))))
    
        #Abweichung h2 
        opth2<-      c(h21, h22, h23, h24, h21c, h22c, h23c, h24c, eh2, mean(h21, h22, h23, h24), mean(h21c, h22c, h23c, h24c))
        abwh2<-      abs(opth2-h2_true)
        ih2<-        which.min(abwh2)
    
        #Abweichung r
        optr<-      c(r_estim, er)
        abwr<-      abs(optr-0.38)
        ir<-        which.min(abwr)
    
        #Abweichung sigHYB
        optsHYB<-      c(sigHYB_estim, sigHYB_estimc, esigHYB)
        abwsHYB<-      abs(optsHYB-0.4^2)
        isHYB<-        which.min(abwsHYB)
    
        #Abweichung sigSNP
        optsSNP<-     c(sigSNP_estim, sigSNP_estimc, esigSNP)
        abwsSNP<-     abs(optsSNP-0.015^2)
        isSNP<-        which.min(abwsSNP)
    
        #Create new Inputvector
        xe[(nf+na+nSNP+1)]<- optr[ir] #0.51      #0.61   #optr[ir]   
        xe[(nf+na+nSNP+2)]<- sqrt(opth2[ih2])    #sqrt(0.3) #sqrt(opth2[ih2]) #sqrt(eh2)  #sqrt(mean(c(h21, h22, h23)))
        xe[(nf+na+nSNP+3)]<- sqrt(optsHYB[isHYB])#0.99      #5.55   #sqrt(mean(optsHYB))   #er r_estim
        xe[(nf+na+nSNP+4)]<- sqrt(optsSNP[isSNP])#0.0080    #0.0046 #sqrt(mean(optsSNP))   #sqrt(esigSNP)
        xe[1]<-abs(xe[1])
        #-----------------------------------------------
        # Scale xe
        #-----------------------------------------------
        scl<<-meSCALE(xe)
        xe1<-  xe/scl
    
        #-----------------------------------------------
        # Solve with optim
        #-----------------------------------------------  
        lb<-          c(rep(-10, nf), rep(-10, na), rep(-10,nSNP), 0.01, 0.1 , 0.01 , 0.0001)/scl
        ub<-          c(rep( 10, nf), rep( 10, na), rep( 10,nSNP), 0.99, 0.99, 1    , 15)/scl
        sol<-         optim(xe1, fmml, grad, method="L-BFGS-B", lower = lb, upper = ub, control = list(maxit=config$maxit2))
    
        print(paste0(    'Iterationen: ',nrow(param[!is.na(param$Rep)])))
        print(Sys.time())
    
        #Output of optim
        print(sol$value)
        print(sol$counts)
        print(sol$convergence)
        print(sol$message)
    
        xout<-            sol$par*scl
    
        er<-              xout[length(xout)-3]
        eh2<-             xout[length(xout)-2]^2
        esigHYB<-         xout[length(xout)-1]^2
        esigSNP<-         xout[length(xout)]^2
        print(paste0(     'h2=', eh2))
        print(paste0(     'r=', er))
        print(paste0(     'sigma_r='  , esigHYB))
        print(paste0(     'sigma_SNP=', esigSNP))

        print(paste0('Korrelation zur vorherigen Loesung=',cor(xe, xout)))
        kor<-         cor(xe, xout)  
    
        if (PRT==1) {
            #Write output for all iterations
            fwrite(as.data.table(param),   paste0(PATH,'.ebv_iter',suffix), col.names=TRUE, row.names=FALSE,quote=FALSE,sep=' ',dec='.')
        }
    
    }
     
#-----------------------------------------------
# Remove unused rows or param
#-----------------------------------------------
    param<-  na.omit(param)
    
#-----------------------------------------------
# Estim parameter for all iterations
#-----------------------------------------------
    estim<-  NULL
    par<-    as.matrix(param)
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
    EBVcan<-  NULL
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
        fwrite(EBVcan,   paste0(PATH,'.solaniTEST',suffix), col.names=TRUE, row.names=FALSE,quote=FALSE,sep=' ',dec='.') 
        #Write output for all iterations
        fwrite(as.data.table(param),   paste0(PATH,'.solITER',suffix), col.names=TRUE, row.names=FALSE,quote=FALSE,sep=' ',dec='.')    
    }
    
    #Plot correlations for test animals
    if (REP>1) {
        plot_corCanREP(EBVcan, param, TRUE)
    } else {
        plot_corCan(EBVcan, param, TRUE)
    }
    
    #Plot output (EBVs fro training animals and parameters)
    if (REP>1) {
        plot_outREP(param, PATH, suffix)
    } else {
        plot_out(param, PATH, suffix)
    }
    
    #Write output files solani, solsnp, solfix and solpar for best solution
    solpar<-           param[ITB1, c('r', 'h2', 'sigSNP', 'sigHYB')]
    solfix<-           param[ITB1, (1+1):nf]
    solani<-           param[ITB1,(nf+1+1):(nf+na)]
    solaniTest<-       EBVcan[,c(1,ITB1)]
    solsnp<-           param[ITB1,(nf+na+1+1):(nf+na+nSNP)]
    fwrite(solpar,     paste0(PATH,'.solpar',suffix),      col.names=TRUE,  row.names=FALSE,quote=FALSE,sep=' ',dec='.')
    fwrite(solfix,     paste0(PATH,'.solfix',suffix),      col.names=FALSE, row.names=FALSE,quote=FALSE,sep=' ',dec='.')
    fwrite(solani,     paste0(PATH,'.solani',suffix),      col.names=FALSE, row.names=FALSE,quote=FALSE,sep=' ',dec='.')
    fwrite(solaniTest, paste0(PATH,'.solani_test',suffix), col.names=FALSE, row.names=FALSE,quote=FALSE,sep=' ',dec='.')
    fwrite(solsnp,     paste0(PATH,'.solsnp',suffix),      col.names=FALSE, row.names=FALSE,quote=FALSE,sep=' ',dec='.')    
    
    print(Sys.time())
    warnings()
q('no')

