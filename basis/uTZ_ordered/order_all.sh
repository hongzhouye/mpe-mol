#!/bin/bash

run_path=$PWD
raw_path=$run_path/../uTZ_qchem

cd $raw_path
element_list=`ls *.txt | sed "s:.txt::g"`
cd $run_path

for I in $element_list
do
	python put_order.py $raw_path/${I}.txt > $I.txt
done
