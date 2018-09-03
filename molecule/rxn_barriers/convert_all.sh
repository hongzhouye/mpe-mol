#!/bin/bash

run_path=$PWD
raw_path=$run_path/../rxn_barriers_raw

cd $raw_path
mols=`ls *.xyz | sed "s:.xyz::g"`
cd $run_path

for mol in $mols
do
	echo "Running for $mol..."
	python xyz2molecule.py $raw_path/$mol.xyz > $mol.in
done
