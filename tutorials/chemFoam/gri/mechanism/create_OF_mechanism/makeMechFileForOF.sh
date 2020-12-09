#!/bin/bash

#Grep the mechanism file of the pyjac in order to retain the right species order


mechPath=../pyJac/out/mechanism.h

#Determine the species within a pyjac output file
grep -A100000 'Species Indexes' $mechPath  > tmp.txt
grep -B100000 '*/' tmp.txt > tmp2.txt
sed '1d;$d' tmp2.txt > tmp3.txt
awk '{ print $2}' tmp3.txt > tmp4.txt
rm -f tmp.txt tmp2.* tmp3.* 

NS=$(cat tmp4.txt | wc -l)


echo "species" >> tmp
echo "$NS" >> tmp
echo "(" >> tmp

cat tmp tmp4.txt > foam_mech_files/chem.foam

rm -f tmp*

printf ");" >> foam_mech_files/chem.foam







