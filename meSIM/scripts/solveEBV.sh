#!/bin/bash
cd ebv

#Delete old files
file="*.log"
if [ -f $file ] ; then
    rm $file
fi

file="Solani.ebv"
if [ -f $file ] ; then
    rm $file
fi

file="Solfix.ebv"
if [ -f $file ] ; then
    rm $file
fi

#Copy mix DIR file and stop file
cp ../scripts/mix_ebv.in  .
cp ../scripts/zws.stop    .

d=$(wc -l 'ebv.dat')
p=$(wc -l 'ebv.ped')
printf "\n%15s" "Phaenotypen:";printf "%8s" $d;printf "\n%15s"  "Pedigree:";printf "%8s" $p;printf "\n"

#Start evaluation
mix99i <mix_ebv.in >mix99i.log; exs=$?
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
mv Solani Solani.ebv
mv Solfix Solfix.ebv

if [ $? -eq 0 ]; then
    printf "\n%30s" "!!!!>>BLUP evaluation done<<!!!!";printf "\n"
else
    printf "\n%30s" "!!!!>>BLUP evaluation failed<<!!!!";printf "\n"
fi



