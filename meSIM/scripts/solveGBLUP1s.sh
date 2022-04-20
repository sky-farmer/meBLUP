#!/bin/bash

cd gblup1s

#Delete old files
file="*.log"
if [ -f $file ] ; then
    rm $file
fi

file="Solani.blup"
if [ -f $file ] ; then
    rm $file
fi

file="Solfix.blup"
if [ -f $file ] ; then
    rm $file
fi

#Copy mix DIR file and stop file
cp ../scripts/mix_1step.in .
cp ../scripts/zws1s.stop   zws.stop


d=$(wc -l 'gblup1s.dat')
p=$(wc -l 'gblup1s.ped')
g=$(wc -l 'animals_geno_sub_all.txt')
printf "\n%15s" "Phaenotypen:";printf "%8s" $d;printf "\n%15s"  "Pedigree:";printf "%8s" $p;printf "\n%15s" "Genotypen:";printf "%8s" $g;printf "\n"


#################################
# G Matrix
###################################
MKL_NUM_THREADS=20 hginv_lapack_para  -m PvR1 -a ../bfq.txt -add 0.001 -z 1E-6 geno_sub_all.txt -lower  Gi.bin > hginv.log.tmp ;exs=$?
if [ $exs -eq 0 ];then printf  "\n%30s" "!!!!>>hginv G setup (full) erfolgreich<<!!!!\n"
                  else printf  "\n%30s" "!!!!>>hginv G setup (full) abgestuerzt<<!!!!\n"
fi

#Start evaluation
mix99i <mix_1step.in >mix99i.log; exs=$?
    if [ $exs -eq 1 ]; then
        printf "\n%30s" "!!!!>>Mix99i  failed<<!!!!";printf "\n"
    else
        printf "\n%30s" "!!!!>>Mix99i  done<<!!!!";printf "\n"
    fi
mix99s <zws.stop >mix99s.log; exs=$?
    if [ $exs -eq 1 ]; then
        printf "\n%30s" "!!!!>>Mix99s  failed<<!!!!";printf "\n"
    else
        printf "\n%30s" "!!!!>>Mix99s  done<<!!!!";printf "\n"
    fi
mv Solani Solani.gblup1s
mv Solfix Solfix.gblup1s

if [ $? -eq 0 ]; then
    printf "\n%30s" "!!!!>>BLUP evaluation done<<!!!!";printf "\n"
else
    printf "\n%30s" "!!!!>>BLUP evaluation failed<<!!!!";printf "\n"
fi


