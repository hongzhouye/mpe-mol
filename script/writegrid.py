from sys import argv
import re
import os
from re_arrange_ao_grid import ao_grid_xform

# Read number of atoms
def atom_number(name):
    f=open(name,'r')
    natom=0; read=0
    while True:
      line=f.readline()
      if "$end" in line:
        break
      if "$molecule" in line:
        f.readline()
        f.readline()
        read=1
      if read==1:
        natom+=1 
    return natom

# Read AO orbital values and gradients on the grid
def ao_grid(name):
    f=open(name,'r')
    o=open('ao_grid.txt','w')
    natom=atom_number(name)
    read=0; atomcount=0
    while True:
      line=f.readline()
      if atomcount>=natom:
        break
      if 'Density 1 and gradient:' in line:
        read=0; atomcount+=1
      if 'Significant basis function' in line:
        read+=1
        f.readline()
        line=f.readline()
      if read>0:
        tmp=re.findall(r"[-+]?\d{0,7}\.\d{8}",line)
        o.write(str(tmp[0])+' '+str(tmp[1])+' '+str(tmp[2])+' '+str(tmp[3])+'\n')

    # @@HY close opened files
    f.close()
    o.close()

name = argv[1]

# root directory where this script is stored
rd = os.path.dirname(os.path.realpath(__file__))
# current working directory
cwd = os.getcwd()

# @@HY: simply call a shell script to get grid dimension and grid itself
# @@HY: note that the grid is now in 4*np layout in order to be processed by C
os.system("bash {:s}/grep_grid_dim.sh {:s}/{:s}_grid.out".\
    format(rd, cwd, name))
ao_grid('{:s}/{:s}_grid.out'.format(cwd, name))
# now the grid in "ao_grid.txt" is in the q-chem style
# which is silly grouped by atoms
# I wrote a separate function to make it suitable for our calculations
ao_grid_xform("{:s}/grid_dimension.txt".format(cwd), \
    "{:s}/ao_grid.txt".format(cwd), \
    "{:s}/ao_grid_new.txt".format(cwd))
os.system("mv {:s}/ao_grid_new.txt {:s}/ao_grid.txt".format(cwd, cwd))
