##########################################
#R-Skript for adding Heard-Year-Effect to phenotypes of simulated data from alphaSimR
#
# USAGE:
#   make_mix99(generation, H2, pop, founder, ZWS)
#       generation....TRUE (erste ZE) or FALSE (>1. ZE)
#       H2............heritability (e.g. 0.4) 
#       pop...........current population (e.g. SELblup)
#       ZWS...........ZWS-type (e.g. blup, gblup1s, gblup2s)
#
# Himmelbauer, 15.6.2021
###########################################

make_mix99<- function (first, H2eval, pop, ZWS) {

#first='FALSE'
#H2=0.3
#pop='SELebv'
#ZWS='ebv'

############################
#Input
############################
ONLY_DAMS_PHENOTYPE=TRUE
#------------------------
#population parameters
PWD<-         paste0(pop, '/')      
H2<-          H2eval   

#------------------------
PWD_ZWS<-     paste0(ZWS,'/')    
DAT<-         ZWS

if (ZWS=='sebv')     { DAT<-  'ebv' }
if (ZWS=='sgblup1s') { DAT<-  'sgblup1s' }
   
#------------------------------------------------------------
#Read Data for current generation
#------------------------------------------------------------
#read Pedigree and Info
ped<-         fread(paste0(PWD,'ped.txt'))
#------------------------
#read true gv
tbv<-         fread(paste0(PWD,'tbv.txt'))
names(tbv)<-  c('ID', 'TBV')
#------------------------
#read pheno
pheno<-       fread(paste0(PWD,'pheno.txt'))
names(pheno)<-c('ID', 'fix', 'trait')

if (first) {

#######################################################################################################
#erste ZE
#######################################################################################################
pedsel<-    ped[ped$ID%in%pheno$ID]; setkey(pedsel,'ID'); setkey(pheno,'ID'); setkey(tbv,'ID')

#--------------------------------------------------------------
#calculate pheno_new
#--------------------------------------------------------------
pheno$trait_new<- pheno$trait
if (ONLY_DAMS_PHENOTYPE==TRUE) {pheno$trait_new[pedsel$SEX=='M']<-NA}

#--------------------------------------------------------------
#Merge data
#--------------------------------------------------------------
data_tmp<-		  data.table(pedsel,phen=pheno$trait_new,phen_orig=pheno$trait, fix=pheno$fix, tbv=tbv$TBV,final_EBV=NA,YrHrd=0,Herdeff=0,Yeareff=0)
fwrite(data_tmp,  paste0(PWD, 'data_tmp.txt'),quote=F,sep=' ',col.names=T,na='NA')

#------------------------------------------------------------
#create files for Mix99-ZWS
#------------------------------------------------------------

#------------------------
#pedigree
fwrite(ped,      paste0(PWD_ZWS,DAT,'.ped') ,quote=F,sep=' ',col.names=F)
fwrite(ped,      paste0(PWD    ,'ped_o.txt'),quote=F,sep=' ',col.names=T)

#------------------------
#phenotypes
last<-            ped[ped$ZE==max(ped$ZE)]
phenon<-          data.table(ID=data_tmp$ID, fix=data_tmp$fix,  pheno=data_tmp$phen)
phenon<-	      phenon[!is.na(phenon$pheno) & !phenon$ID%in%last$ID,]
fwrite(phenon,    paste0(PWD_ZWS,DAT, '.dat'),quote=F,sep=' ',col.names=F)

#------------------------
#mix.var file based on heritabilities
sink(paste0(PWD_ZWS,DAT,'.var'))
	cat(paste(1,' ',1,' ',1,' ',H2,sep=''),'\n')
	cat(paste(2,' ',1,' ',1,' ',1-H2,sep=''),'\n')
sink()

}

if (!first) {
##################################################################################################
#Generation>0
##################################################################################################
lastZe<-       ped[ped$ZE==max(ped$ZE)]
tbv<-          tbv[tbv$ID%in%lastZe$ID]
pheno<-        pheno[pheno$ID%in%lastZe$ID]

#------------------------------------------------------------
#Read old Data
#------------------------------------------------------------
data_o<-	      fread(paste0(PWD, 'data_tmp.txt'), header=T,na.strings='NA')
ped_o<-           fread(paste0(PWD, 'ped_o.txt'), sep=' ',header=T)

#--------------------------------------------------------------
#calculate pheno_new
#--------------------------------------------------------------
pheno$trait_new<- pheno$trait
if (ONLY_DAMS_PHENOTYPE==TRUE) {pheno$trait_new[lastZe$SEX=='M']<-NA}

#--------------------------------------------------------------
#Merge data
#--------------------------------------------------------------
data_tmp<-		  data.table(lastZe,phen=pheno$trait_new,phen_orig=pheno$trait, fix=pheno$fix, tbv=tbv$TBV,final_EBV=NA,YrHrd=0,Herdeff=0,Yeareff=0)

#------------------------
#merge with old data and data.tmp file is written containing all information on pedigree, sex, generation, true genotype	
data_tmp<-	      rbind(data_o,data_tmp)
fwrite(data_tmp,  paste0(PWD, 'data_tmp.txt'),quote=F,sep=' ',col.names=T,na='NA')

#------------------------------------------------------------
#create new data for Mix99-ZWS
#------------------------------------------------------------

#------------------------
#pedigree
new<-             ped[!ped$ID%in%ped_o]
pedn<-		      unique(rbind(ped_o,new)); setkey(pedn, 'ID')
fwrite(pedn,      paste0(PWD_ZWS,DAT, '.ped') ,quote=F,sep=' ',col.names=F)
fwrite(pedn,      paste0(PWD    , 'ped_o.txt'),quote=F,sep=' ',col.names=T)

#------------------------
#phenotypes
last<-            ped[ped$ZE==max(ped$ZE)]
phenon<-          data.table(ID=data_tmp$ID, fix=data_tmp$fix, pheno=data_tmp$phen)
phenon<-	      phenon[!is.na(phenon$pheno) & !phenon$ID%in%last$ID,]
fwrite(phenon,    paste0(PWD_ZWS,DAT, '.dat'),quote=F,sep=' ',col.names=F)
}

}





