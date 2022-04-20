####################################################
# Simulation mit AlphaSimR
#
#
# Himmelbauer, 16.06.2021 
####################################################

print(Sys.time())

#---------------------------------------------------------------------
#Package check
#---------------------------------------------------------------------
packages = c("data.table", "AlphaSimR")   # packages needed in the Script

# Load or install&load all packeges
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

#---------------------------------------------------------------------
#Load functions
#---------------------------------------------------------------------
source('scripts/make_mix99.r')

#---------------------------------------------------------------------
#Define bash-scripts
#---------------------------------------------------------------------
solveEBV<-    'scripts/solveEBV.sh'      #konv Schaetzung
solveGBLUP1s<-'scripts/solveGBLUP1s.sh'  #single Step Schaetzung
solveEBV2<-   'scripts/solveEBV2.sh'     #konv Schaetzung 2

##########################
# Global parameters
##########################
createFounder=FALSE  #Create founder population or load it from data
nfounder = 2000  #Founder animals in historic population in generation 0
nQtl     = 30    #Number of QTL per chromosome
nSnp     = 1660  #Number of SNP per chromosome
nCHR     = 30    #Number of chromosomes
nPopBn   = 150   #effective populationsize
BP       = 1e+08 #Number of base pairs per chromosome
genL     = 1     #Genetic length in Morgan per chromosome
#Trait
initMeanG   = 0    #Initial mean
initVarG    = 1    #Initial variances Genetic
H2          = 0.5  #Heritability

#SelEBV
nSires      =50
nDams       =1000
nGEN        =11
nSires_old  =5
nDams_old   =700

##########################
# Generate initial haplotypes
##########################
if (createFounder==TRUE) {
FOUNDERPOP = runMacs2(  nfounder,
                        nChr     = nCHR,
                        segSites = nQtl+nSnp,
                        Ne       = nPopBn,
                        bp       = BP,
                        genLen   = genL,
                        mutRate  = 2.5e-08,   #default 
                        inbred   = FALSE,
                        split    = NULL,
                        ploidy   = 2L,
                        nThreads = 20
                      )
#-----------------------------------
# Save FOUNDERPOP to file
saveRDS(FOUNDERPOP, "FOUNDERPOP.rds")
#-----------------------------------
} else {
FOUNDERPOP <- readRDS("FOUNDERPOP.rds")
}

#-----------------------------------
#Set simulation parameters
SP = SimParam$new(FOUNDERPOP)
SP$restrSegSites(nQtl,nSnp)
if(nSnp>0){
  SP$addSnpChip(nSnp)
}

SP$addTraitA(nQtl,mean=initMeanG,var=initVarG)
SP$setVarE(h2=H2)       
SP$setSexes("yes_sys")
SP$setTrackPed(TRUE)
SP$setTrackRec(TRUE)

#---------------------------------------------------------------------
# Founder population
#---------------------------------------------------------------------
hp= newPop(FOUNDERPOP)
hp= setPheno(hp, h2=H2_1)

m<-        selectInd(pop=hp, nSires,use="rand",sex="M",returnPop=TRUE,simParam=SP)  #all selected males  rep((ze-1),nInd(old))
f<-        selectInd(pop=hp, nDams,use="rand",sex="F",returnPop=TRUE,simParam=SP)  #all selected females     
m@iid<-    c(1:nSires)
m@id<-     as.character(c(1:nSires))
f@iid<-    c((nSires+1):(nSires+nDams))
f@id<-     as.character(c((nSires+1):(nSires+nDams)))

#Founder population
hp<-       c(m,f)

#Save allele frequencies
geno<- pullSnpGeno(hp)
sum<-  colSums(geno)
nanim<-nrow(geno)
fq<-   data.table(SNP<-1:length(sum), f=sum/(2*nanim))
fwrite(fq, "bfq.txt", col.names=FALSE, row.names=FALSE, quote=FALSE, sep=' ', dec='.', append=FALSE)


##########################
#Populations
##########################

#-------------------------------------------------
#SELebv
#-------------------------------------------------
old<-              hp
old@misc$gen<-     rep(0,nInd(old))  #Add gernation number

CHK<-              TRUE
CHK2<-             FALSE

for(ze in 1:nGEN){
    new = randCross(pop=old,
                    nCrosses=1000,         #CHANGE to 30000
                    nProgeny = 1,
                    simParam = SP,
                    balance = TRUE)
    
    #Add generation number
    new@misc<-     rep(list((ze)),nInd(new)) 
    
    #Write output
    ped<-          data.table(ID=new@id, SIRE=new@father, DAM=new@mother, SEX=new@sex, ZE=unlist(new@misc))
    pheno<-        data.table(ID=new@id, fixEff=new@fixEff, pheno=new@pheno)
    tbv<-          data.table(ID=new@id, tbv=new@gv)
    if (ze==1) {
        ped<-      rbind(data.table(ID=old@id, SIRE=old@father, DAM=old@mother, SEX=old@sex, ZE=unlist(old@misc)), ped)
        pheno<-    rbind(data.table(ID=old@id, fixEff=old@fixEff, pheno=old@pheno), pheno)
        tbv<-      rbind(data.table(ID=old@id, tbv=old@gv), tbv)
    }       
    fwrite(ped,   "SELebv/ped.txt",   col.names=CHK, row.names=FALSE,quote=FALSE,sep=' ',dec='.', append=CHK2)
    fwrite(pheno, "SELebv/pheno.txt", col.names=CHK, row.names=FALSE,quote=FALSE,sep=' ',dec='.', append=CHK2)
    fwrite(tbv,   "SELebv/tbv.txt",   col.names=CHK, row.names=FALSE,quote=FALSE,sep=' ',dec='.', append=CHK2)
 
    #Save genotype
    if (ze >= 3) {
        genoall<-     data.table(ID=new@id, as.data.table(pullSnpGeno(new)))
        anim_genoall<-data.table(ID=new@id)
        if (ze==3) { CHK3=FALSE } else { CHK3=TRUE }
        fwrite(genoall,      "SELebv/geno_all.txt",       col.names=FALSE, row.names=FALSE,quote=FALSE,sep=' ',dec='.', append=CHK3)
        fwrite(anim_genoall, "SELebv/anims_geno_all.txt", col.names=FALSE, row.names=FALSE,quote=FALSE,sep=' ',dec='.', append=CHK3)
    }
    
    #Create files for estimation with mix99
    make_mix99(CHK, H2_1, 'SELebv', 'ebv')
    
    #Esimate EBV with Mix99 BLUP
    system(solveEBV,wait=T)
    
    #Import and set EBVs
    EBV<-  fread("ebv/Solani.ebv", colClasses=c('character', 'NULL', 'NULL', 'numeric')); setkey(EBV, 'V1')
    EBV_o<-EBV[J(old@id), nomatch=NA]
    EBV_n<-EBV[J(new@id), nomatch=NA]
    old@ebv = as.matrix(EBV_o$V4) #Assign EBVs to population 
    new@ebv = as.matrix(EBV_n$V4) #Assign EBVs to population 
    
    #select individuals from old generation
    male_old<-     selectInd(pop=old, nSires_old, use = "ebv", selectTop = TRUE, sex = "M", returnPop = TRUE, simParam = SP)          #CHANGE to 75
    female_old<-   selectInd(pop=old, nDams_old, use = "ebv", selectTop = TRUE, sex = "F", returnPop = TRUE, simParam = SP)          #CHANGE to 21000
   
    #select individuals from new generation
    male_new<-     selectInd(pop=new,(nSires-nSires_old), use = "ebv", selectTop = TRUE, sex = "M", returnPop = TRUE, simParam = SP)          #CHANGE to 675
    female_new<-   selectInd(pop=new,(nDams-nDams_old),   use = "ebv", selectTop = TRUE, sex = "F", returnPop = TRUE, simParam = SP)          #CHANGE to 9000
    sel<-          c(male_new, female_new)
    sel@misc$gen<- rep(ze,nInd(sel))
        
    old<-          c(male_new, female_new, male_old, female_old)
    print(paste0("Done for Timeunit: ",ze))
    CHK<-          FALSE
    CHK2<-         TRUE
}

#Prepare datafiles for meBLUP

#Read Datafiles
ped<-   fread('SELebv/ped.txt'); setkey(ped, 'ID')
dat<-   fread('SELebv/pheno.txt'); setkey(dat, 'ID')
tbv<-   fread('SELebv/tbv.txt'); setkey(tbv, 'ID')
geno<-  fread('SELebv/geno_all.txt'); setkey(geno, 'V1')
agt<-   fread('SELebv/anims_geno_all.txt'); setkey(agt, 'V1')

#Recode animal IDs
recode<-     data.table(oID=ped$ID, nID=1:nrow(ped))
recode<-     rbind(data.table(oID=0, nID=0), recode)
setkey(recode, 'oID')

#ped
setkey(ped, 'ID')
pedID<-      recode[J(ped$ID), nomatch=0]
pedS<-       recode[J(ped$SIRE), nomatch=0]
pedD<-       recode[J(ped$DAM), nomatch=0]
ped<-        data.table(ID=pedID$nID, SIRE=pedS$nID, DAM=pedD$nID, ZE=ped$ZE, SEX=ped$SEX); setkey(ped, 'ID')

#pheno
setkey(dat, 'ID')
ID<-         recode[J(dat$ID), nomatch=0]
dat$ID<-ID$nID; setkey(dat, 'ID')

#tbv
setkey(tbv, 'ID')
ID<-         recode[J(tbv$ID), nomatch=0]
tbv$ID<-  ID$nID; setkey(tbv, 'ID')

#geno
setkey(geno, 'V1')
ID<-         recode[J(geno$V1), nomatch=0]
geno$V1<-    ID$nID

setkey(agt, 'V1')
ID<-         recode[J(agt$V1), nomatch=0]
agt$V1<-     ID$nID


#Select phenotypes
pedpheno<-    ped[ped$ZE<(6) & ped$ZE!=0 & ped$SEX=='F',]
pheno<-       dat[dat$ID%in%pedpheno$ID,]
phenoME<-     data.table(ID=pheno$ID, pheno=pheno$pheno, fixEff=pheno$fix)

#Write recoded Output to SElebv
fwrite(ped,         "SELebv/ped.txt",           col.names=TRUE, row.names=FALSE,quote=FALSE,sep=' ',dec='.', append=FALSE)
fwrite(dat,         "SELebv/pheno.txt",         col.names=TRUE, row.names=FALSE,quote=FALSE,sep=' ',dec='.', append=FALSE)
fwrite(tbv,         "SELebv/tbv.txt",           col.names=TRUE, row.names=FALSE,quote=FALSE,sep=' ',dec='.', append=FALSE)
fwrite(geno,        "SELebv/geno_all.txt",      col.names=TRUE, row.names=FALSE,quote=FALSE,sep=' ',dec='.', append=FALSE)
fwrite(agt,         "SELebv/anims_geno_all.txt",col.names=TRUE, row.names=FALSE,quote=FALSE,sep=' ',dec='.', append=FALSE)

#Write Output for ssGBLUP
fwrite(ped,         "gblup1s/gblup1s.ped",              col.names=FALSE, row.names=FALSE,quote=FALSE,sep=' ',dec='.', append=FALSE)
fwrite(pheno,       "gblup1s/gblup1s.dat",              col.names=FALSE, row.names=FALSE,quote=FALSE,sep=' ',dec='.', append=FALSE)
fwrite(geno,        "gblup1s/geno_sub_all.txt",         col.names=FALSE, row.names=FALSE,quote=FALSE,sep=' ',dec='.', append=FALSE)
fwrite(agt,         "gblup1s/animals_geno_sub_all.txt", col.names=FALSE, row.names=FALSE,quote=FALSE,sep=' ',dec='.', append=FALSE)
system('cp ebv/ebv.var   gblup1s/gblup1s.var')

#Write Output for meBLUP
ped<-ped[,-5]
fwrite(ped,         "meBLUP/meBLUP.ped",       col.names=FALSE, row.names=FALSE,quote=FALSE,sep=' ',dec='.', append=FALSE)
fwrite(phenoME,     "meBLUP/meBLUP.dat",       col.names=FALSE, row.names=FALSE,quote=FALSE,sep=' ',dec='.', append=FALSE)
fwrite(geno,        "meBLUP/meBLUP.snp",       col.names=FALSE, row.names=FALSE,quote=FALSE,sep=' ',dec='.', append=FALSE)
fwrite(agt,         "meBLUP/animals_geno.txt", col.names=FALSE, row.names=FALSE,quote=FALSE,sep=' ',dec='.', append=FALSE)
system('cp ebv/ebv.var   meBLUP/meBLUP.var')

system(solveGBLUP1s)  #single Step Schaetzung
#system(solveEBV2)     #konv Schaetzung

#Merge Ergebnisse
#ebv2<-    fread('ebv2/Solani.ebv2');      setkey(ebv2,    'V1'); ebv2<-ebv2[ebv2$V1>0]
gblup1s<- fread('gblup1s/Solani.gblup1s');setkey(gblup1s, 'V1'); gblup1s<-gblup1s[gblup1s$V1>0]
tbv<-     fread('SELebv/tbv.txt');        setkey(tbv, 'ID')
ped<-     fread('SELebv/ped.txt');        setkey(ped, 'ID')
bv<-      data.table(id=tbv$ID, gen=ped$ZE ,true=tbv$tbv.V1, sstep=gblup1s$V4)
fwrite(bv,         "meBLUP/meBLUP.bv",     col.names=TRUE, row.names=FALSE,quote=FALSE,sep=' ',dec='.', append=FALSE)

print(Sys.time())
q('no')








