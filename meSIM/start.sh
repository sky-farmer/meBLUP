#!/bin/bash
#echo "start.sh" | at now

#######################################
STARTFILE='scripts/startSimu.r'
######################################

#Create subdirectories
if [ ! -d 'SELebv' ];     then  mkdir SELebv    ; fi
if [ ! -d 'ebv' ];        then  mkdir ebv       ; fi
if [ ! -d 'ebv2' ];       then  mkdir ebv2      ; fi
if [ ! -d 'SELgblup1s' ]; then  mkdir SELgblup1s; fi
if [ ! -d 'gblup1s' ];    then  mkdir gblup1s   ; fi
if [ ! -d 'gblup1s/si' ]; then  mkdir gblup1s/si; fi

#Starte Simulation
R CMD BATCH  $STARTFILE
    
#Loesche Ueberfluessige Dateien	
rm SELgblup1s/geno.txt
rm SELgblup1s/ped_o.txt
rm SELgblup1s/pheno.txt
rm SELebv/ped.txt
rm SELebv/pheno.txt
rm SELebv/tbv.txt
rm SELebv/geno_all.txt
rm SELebv/anims_geno_all.txt
rm SELebv/data_tmp.txt
rm SELebv/ped_o.txt
rm ebv/*
rm gblup1s/Gi.bin
rm gblup1s/geno_sub_all.txt
