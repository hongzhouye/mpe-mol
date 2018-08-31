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

cd $path

# obtain grid
/home/hzye2011/local/bin/shutils/qchem.20180508 ${name}_grid.in ${name}_grid.out

python ${run_path}/writegrid.py $name
