#!/bin/bash

cd ebv2

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
cp ../scripts/mix_ebv2.in .
cp ../scripts/zws1s.stop   zws.stop


d=$(wc -l '../gblup1s/gblup1s.dat')
p=$(wc -l '../gblup1s/gblup1s.ped')
printf "\n%15s" "Phaenotypen:";printf "%8s" $d;printf "\n%15s"  "Pedigree:";printf "%8s" $p;printf "\n%15s" "Genotypen:";printf "%8s" $g;printf "\n"


#Start evaluation
mix99i <mix_ebv2.in >mix99i.log; exs=$?
if [ $exs -eq 1 ]; then
	printf "\n%30s" "!!!!>>Mix99i (conv.)  failed<<!!!!";printf "\n"; exit 1
else
    printf "\n%30s" "!!!!>>Mix99i (conv.)   done<<!!!!";printf "\n"
fi

mix99s <zws.stop >mix99s.log; exs=$?
    if [ $exs -eq 1 ]; then
        printf "\n%30s" "!!!!>>Mix99s (conv.)   failed<<!!!!";printf "\n"; exit 1
    else
        printf "\n%30s" "!!!!>>Mix99s (conv.)   done<<!!!!";printf "\n"
    fi


mv Solani Solani.ebv2
mv Solfix Solfix.ebv2 




