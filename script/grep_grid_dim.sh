#!/bin/bash

# goal1: cat Natom Npoints Nbasis > grid_dimension.txt
# goal2: cat GRIDPOINTS > grid.txt
# usage: ./grep_grid_dim.sh bh_grid.out

if [ $# -lt 1 ]
then
	echo "Usage: Q-Chem_output"
	exit
fi

inp=$1
prefix=`echo $inp | sed "s:.out::g"`

grepQChemStructure.sh $prefix.out 2 "last" > /dev/null
natom=`wc -l ${prefix}.xyz | awk '{print $1-2}'`
rm ${prefix}.xyz

# npoint=`grep -a "^Grid " $prefix.out | wc -l | awk '{print $1}'`
# the statement above only gives total number of grid points
# but since the grid points are grouped by atom number in (silly) q-chem outputs
# we need the grid point for each atom
grep -an "Grid 0:" $prefix.out > _hy_tmp1
grep -an "Significant basis function  1 " $prefix.out > _hy_tmp2
for I in `seq 1 1 $natom`
do
	nstart=`awk "NR==$I{print}" _hy_tmp1 | awk -F ":" '{print $1}'`
	nend=`awk "NR==$I{print}" _hy_tmp2 | awk -F ":" '{print $1}'`
	npoints[$I]=`echo $nstart $nend | awk '{print $2-$1}'`	
done
rm _hy_tmp1 _hy_tmp2

nbasis=`grepQChemBasOcc.sh $prefix.out | grep -a "Nbasis" | awk '{print $3}'`

echo $natom ${npoints[*]} $nbasis > grid_dimension.txt

# grep grid points info and transpose them in order to match C I/O style
# the reason I did this was that I would like to ask C to do all grid-related operations
# simply because I DON'T KNOW HOW TO CODE IN FORTRAN!!!
# the script for transposition was adapted from a stack overflow answer:
# https://stackoverflow.com/questions/1729824/an-efficient-way-to-transpose-a-file-in-bash
grep -a "^Grid" $prefix.out | awk '{print $3" "$4" "$5" "$6}' | \
awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' > grid.txt
