################################################################################
#################################meREAD.r#######################################
################################################################################
# R-Script/function to read Data for meBLUP
#
#-------------------------------------------
#USAGE:
# meRead('DATASET','PRT')
#   DATASET....path to input files + name of the dataset
#   PRT........printlevel (=1)
#
#--------------------------------------------
#INPUT:
#   DATASET.ped........pedigree file [ID SIRE DAM GENERATION] in addition if available BVs [tbv  ssBLUP meBLUP] (missing parent=0)  
#   DATASET.dat........phenotypic data file [ID PHENO FIX(=1)]
#   DATASET.snp........genotype data file [ID GENOTYPES(=012 coded without blanks in beetween)]
#
#OUTPUT:
#  y..........data vector of model equation 
#  A..........data matrix  of model equation 
#       .S.........large, dense matrix of SNPs in char format (A.S-49 turns characters into integers -1,0,1)
#       .X.........dense fixed effect matrix with few columns
#       .PA........sparse parent average matrix for training animals (for defining P)
#       .PAcan.....sparse parent average matrix for test animals
#       .as........list of animals with genotype  (for defining P)
#       .ap........list of animals with phenotype (for defining Z)
#       .ac........list of animal classes (for defining C_H)
#       .can.......list of candidate animals (test animals)
#       .cal.......list of calibration animals (training animals)
#  classDef...class definition, [ac nc gen Ts Ss Ds]
#               (0: not genotyped, 1: genotyped, 2: unknown parent)
# allBV.......[TBV EBVss EBVme] breeding values for comparison
#           TBV......true breeding values
#           EBVss....estimated breeding values from ssBLUP
#           EBVme....estimated breeding values from meBLUP
#
# Himmelbauer, 14.05.2021 (nach Neumaier 2021, meRead.m) | last update 19.04.2022 
##################################################################################

library(data.table)
library(Matrix)


meRead=function(dataset,prt){
    DATASET=dataset
    if (prt) {
    print(DATASET)
    }
    
    A<-  list()
    
    
    ##########################################
    # Read and process pedigree information 
    #     A.PA....contains sparse parent average matrix (for defining P)
    ##########################################
    
    ped<-        fread(paste0(DATASET,'.ped'), fill=TRUE, header=FALSE)
    na<-         nrow(ped)
    if (!all(((1:na)-ped$V1)==0)) {     #animals are assumed to be labelled consecutively by 1:na
        stop('labels of animals must be 1:na')
    }
    rows<-       data.table(c(ped$V1,ped$V1))
    cols<-       data.table(c(ped$V2,ped$V3)) 
    rows<-       rows$V1[cols$V1!=0]
    cols<-       cols$V1[cols$V1!=0]
    nr<-         length(rows)
    vals<-       0.5+rep(0,nr)
    A$PA<-       sparseMatrix(i=rows,j=cols,x=vals,dims=list(na,na))
    gen<-        ped$V4     #generation number
    if (ncol(ped)==5) { TBV<-    ped$V5;  allBV<-  data.table(TBV)}                 #true breeding value
    if (ncol(ped)==6) { EBVsst<- ped$V6;  allBV<-  data.table(TBV,EBVsst)}          #estimated breeding values from sstBLUP
    if (ncol(ped)==7) { EBVme<-  ped$V7;  allBV<-  data.table(TBV,EBVsst,EBVme) }   #estimated brreding values from meBLUP
    
    if (prt) {
        print('meRead: pedigree assigned!')
    }
    
    
    ##########################################
    # Read and process phenotype information 
    #     A.ap....list of animals with phenotypes
    #     A.X.....fixed effect (pheno mean) ???
    #     y(j)....phenotype of animals ap(j)
    ##########################################
    
    dat<-        fread(paste0(DATASET,'.dat'), fill=TRUE, header=FALSE)
    genID<-      ped[ped$V1 %in% dat$V1, c("V4")]
    pheno<-      dat$V3     #phenotype
    

 
    #check if last generation contains relevant information
    last<-        max(gen)
    ylast<-       dat$V3[genID$V4==last]
    if (length(ylast>0)) {
    MMMylast<-    c(min(ylast), mean(ylast), max(ylast))
    if (!all(MMMylast==0)) {
        print(MMMylast)
        stop('Last generation contains relevant information')
    }
    }
    #ignore entries for last generation
    ap<-         dat$V1[genID$V4<last & dat$V3!=MISS]    #phenotypes except in last generation
    nap<-        length(ap)             #number of animals with phenotypes
    A$ap<-       ap                     #list of animals with phenotypes
    A$X<-        matrix(-999,nrow=nap, ncol=(ncol(dat)-2))
    for (i in 2:(ncol(dat)-1)) {
    c<-          i-1
    sel<-        paste0("V",i)
    A$X[,c]<-    dat[dat$V1%in%ap,get(sel)]         #fixed effect: pheno mean
    }
    y<-          dat$V3[dat$V1%in%ap]   #phenotypes of animals in ap
    
    if (prt) {
        print('meRead: phenotypes assigned!')
    }
    
    ##########################################
    # Candidates are not used for optimization
    #  Candidates are all animals without phenotype and not appearing in pedigree as sire or dam
    ##########################################
    cands<-   ped[!ped$V1%in%ped$V2 & !ped$V1%in%ped$V3 & ped$V4==max(ped$V4)]
    cands<-   cands[!cands$V1%in%A$ap]
    A$can<-   cands$V1
    A$cal<-   ped$V1[!ped$V1%in%cands$V1]
    naALL<-   na
    na<-      length(A$cal)
    
    if (!all(((1:na)-A$cal)==0)) {     #animals are assumed to be labelled consecutively by 1:na
        stop('Animals in calibration must be labelled consecutively by 1:na')
    }
    
    A$PAcan<- A$PA[A$can,A$cal]
    A$PA<-    A$PA[A$cal,A$cal]
    
    gen<-     ped$V4[ped$V1%in%A$cal]
    
    
    ##########################################
    # Read and process genotype information 
    #     A.S.....contains SNPs of animals as in char format
    #     A.as....contains list of genotyped animals
    ##########################################
    
    GENO<-       fread(paste0(DATASET,'.snp'))
    asnp<-       GENO$V1
    
    A$as<-       asnp[asnp%in%A$cal]
    A$ascan<-    asnp[asnp%in%A$can]
    
    A$S<-        as.matrix(GENO[GENO$V1%in%A$cal,-1])
    A$Scan<-     as.matrix(GENO[GENO$V1%in%A$can,-1]) 
    
    if (prt) {
        print('meRead: genotypes assigned!')
    }   
    
    if (NSNP==0) {
    A$as<-       NULL
    A$ascan<-    NULL
    
    A$S<-        NULL
    A$Scan<-     NULL
    
    if (prt) {
        print('meRead: no genotypes used!')
    }
    }
    
    ##########################################
    # Create animal classes 
    #     A.S.....contains SNPs of animals as in char format
    #     A.as....contains list of genotyped animals
    #     A.ac(T)...class of animal T in generation gen with parents s and d
    #               gen....generation of T
    #               ts=1...if T genotyped, 0 otherwise
    #               vs=2...if V unknown, 1 if V genotyped, 0 otherwise
    #               ms=2...if M unknown, 1 if M genotyped, 0 otherwise 
    ##########################################
    
    Ts<-         rep(0,na); Ts[A$as]<-1
    Ss<-         rep(0,na)
    Ds<-         rep(0,na)
    for (T in (1:na)) {
        s<-      ped$V2[ped$V1==T]
        if (s==0) {
            Ss[T]<-2
        } else {
            Ss[T]<-Ts[s]
        }
        d<-      ped$V3[ped$V1==T]
        if (d==0) {
            Ds[T]<-2
        } else {
            Ds[T]<-Ts[d]
        }  
    }
    gTSD<-       data.table(gen, Ts, Ss, Ds)
    acode<-      27*gen+9*Ts+3*Ss+Ds+1
    ac<-         0
    Cac<-        NULL
    Cnc<-        NULL
    CgTSD<-      NULL
    aci<-        rep(-999,na)
    for (i in (min(acode):max(acode))) {
        idx<-    which(acode==i)
        nc<-     length(idx)
        if (nc==0) {
            next
        }
        ac<-     ac+1
        Cac<-    rbind(Cac,ac)
        Cnc<-    rbind(Cnc,nc)
        CgTSD<-  rbind(CgTSD, gTSD[idx[1],])
        aci[idx]<-ac
           
    }
    A[["ac"]]<-aci
    classDef<-   data.table(Cac, Cnc, CgTSD)
    names(classDef)<-c('class', 'n', names(CgTSD))
    
    if (prt) {
        print('meRead: animal classes assigned!')
        print(classDef)
        print('0: not gt, 1: gt, 2: unknown parent')
    }   

    
    ##########################################
    # Data consistency checks
    ##########################################

    if (!all(((1:naALL)-ped$V1)==0)) {     #animals are assumed to be labelled consecutively by 1:na
        stop('Animals in pedigree must be sorted')
    }

    if (!all((asnp[2:length(asnp)]-asnp[1:(length(asnp)-1)])>0)) {     #animals are assumed to be labelled consecutively by 1:na
        stop('Animals in genotype file must be sorted')
    }
    
    if (prt) {
        print('meRead: data checks completed')
    }

    #OUTPUT
    if (ncol(ped)>4) {
    return(list(y,A, classDef, allBV))
    } else {
    return(list(y,A, classDef))
    }
}














