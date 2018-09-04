#!/bin/bash

if [ $# != 1 ]
then
	echo "Usage: xxx.out"
	exit
fi

inp=$1

echo "Determining bases dimensions and # of electrons..."
Nbas=`grep -a "shells and" $inp | head -n1 | awk '{print $6}'`
Naux=`grep -a "shells and" $inp | tail -n1 | awk '{print $6}'`
Nele=`grep -a "electrons" $inp | head -n1 | awk '{print $3}'`
echo "$Nbas $Nele $Naux" > dimensions.txt

# grep the block of Ekin and Enuc matrices for python to read
echo "Determining line #'s for Ekin and Enuc matrices..."
Nline_kin=`grep -an "Kinetic Energy Matrix" $inp | awk -F ":" 'NR==1{print $1}'`
Nline_nuc=`grep -an "Nuclear Attraction Matrix" $inp | awk -F ":" 'NR==1{print $1}'`
Nline_dip=`grep -an "Multipole Matrix (1,0,0)" $inp | awk -F ":" 'NR==1{print $1}'`

echo "Dumping Ekin matrix to file..."
awk "NR>=${Nline_kin}&&NR<${Nline_nuc}{print}" $inp > _hy_tmp_kin

echo "Dumping Enuc matrix to file..."
awk "NR>=${Nline_nuc}&&NR<${Nline_dip}{print}" $inp > _hy_tmp_nuc
