#!/bin/bash

#first argument of script is molecule name
#run ./mpeinput.sh name
#adjust scratch directory name

if [ $# != 2 ]
then
	echo "Usage: path prefix"
	exit
fi

run_path=$PWD

path=$1
name=$2

scratchdir=/scratch/hzye2011/$name
qchem=/home/pdesilva/software/qchem.20150304/qchem
#qchem=qchem.latest

cd $path

rm $scratchdir/*
$qchem -save ${name}.in ${name}.out $name

#scratch files used by Qchem, in binary format:
# 22.0 -> 2-electron repulsion integrals in MO basis, only non-zero integrals!
# 809.0 -> MO indices of non-zero integrals in file 22.0
# 53.0 -> MO coefficients followed by MO energies
# 362.0 -> inverse square root of repulsion integrals between auxiliary (RI) basis functions
# 363.0 -> 3-center repulsion integrals between MO pairs and auxiliary (RI) basis functions

mv $scratchdir/22.0 ./
mv $scratchdir/809.0 ./
mv $scratchdir/802.0 ./
mv $scratchdir/53.0 ./
mv $scratchdir/320.0 ./
mv $scratchdir/362.0 ./
mv $scratchdir/363.0 ./

#convert binary files to decimal or integer format

hexdump -v -e '"%.16e "' -e '"\n"' 22.0 > 22.0.txt
hexdump -v -e '4/4 "%d "' -e '"\n"' 809.0 > 809.0.txt
hexdump -v -e '"%.16e "' -e '"\n"' 802.0 > 802.0.txt
hexdump -v -e '"%.16e "' -e '"\n"' 53.0 > 53.0.txt
hexdump -v -e '"%.16e "' -e '"\n"' 320.0 > 320.0.txt
hexdump -v -e '"%.16e "' -e '"\n"' 362.0 > 362.0.txt
hexdump -v -e '"%.16e "' -e '"\n"' 363.0 > 363.0.txt
rm *.0

/home/hzye2011/local/bin/shutils/qchem.20180508 ${name}_grid.in ${name}_grid.out

python ${run_path}/writedata.py $name
python ${run_path}/writegrid.py $name
