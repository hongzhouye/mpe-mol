import numpy as np
import math
from scipy import linalg
from sys import argv
import re

# from a Qchem output file, read the number of orbitals and auxiliary basis functions
# the number of auxiliary basis functions is found by searching for 'rimp2' in output -> change this if you do not use an rimp2 auxiliary basis
def read_number_mos(fname):

    f = open(fname)
    line = f.readline()

    while(not 'basis functions' in line):
        line = f.readline()

    nbf = int(line.rstrip().split()[5])

    # while(not 'basis set is non-standard' in line):
    #     line = f.readline()
    # @@HY
    while(not 'Thank you ' in line):
        line = f.readline()

    line = f.readline()
    while(not 'basis set is non-standard' in line):
        line = f.readline()

    line = f.readline()
    while(not 'basis set is non-standard' in line):
        line = f.readline()

    line = f.readline()
    nribf = int(line.rstrip().split()[5])

    f.close()
    return nbf, nribf

# from a Qchem output file, read the number of occupied orbitals = number of electron pairs
def read_number_nocc(fname):

    f = open(fname)
    line = f.readline()

    while(not 'There are' in line):

        line = f.readline()

    nalpha = int(line.rstrip().split()[2])
    nbeta = int(line.rstrip().split()[5])
    if (nalpha == nbeta):
        nocc = nalpha
    else:
     print "System not spin-compensated!!!"

    f.close()
    return nocc

def read_atom_number(fname):
    f = open(fname)
    line = f.readline()

    while(not '0 1' in line):

        line = f.readline()

    line = f.readline();i=0
    while line.rstrip().split()[0]!='$end':
        i=i+1
        line = f.readline()
    return i

# from a Qchem output file, read the number of atoms
def read_atom_num(fname):
    f = open(fname)
    num_atom=0
    line = f.readline()

    while(not '0 1' in line):

        line = f.readline()

    while(not '$end' in line):
        line = f.readline()
        num_atom+=1

    f.close()
    return num_atom-1

# from a Qchem output file, read atom coordinates
def read_atom_coordinate(fname,num_atom):

    coor = np.zeros((num_atom,3))
    atom = []
    f = open(fname)
    line = f.readline()

    while(not '0 1' in line):

        line = f.readline()

    for i in range(num_atom):
        line = f.readline()
#        print line.rstrip().split()
        atom.append(line.rstrip().split()[0])
        coor[i,0] = float(line.rstrip().split()[1])
        coor[i,1] = float(line.rstrip().split()[2])
        coor[i,2] = float(line.rstrip().split()[3])

    f.close()
    return atom, coor


# Read auxiliary basis set from the output file

def read_aux_basis(fname):

    f=open(fname)

    line = f.readline()
    while(not 'basis set is non-standard' in line):
        line = f.readline()

    line = f.readline()
    while(not 'basis set is non-standard' in line):
        line = f.readline()

    line = f.readline()
    while(not 'basis set is non-standard' in line):
        line = f.readline()

    line = f.readline()
    nshells = int(line.rstrip().split()[2])

    for i in range(7): 
        line = f.readline()
    
    auxbas = []
    for i in range(nshells):
        line = f.readline()
	if (len(line.split())==5):
	    atom = int(line.split()[0])
            auxbas.append(map(int,line.split()[0:3])+map(float,line.split()[3:5]))
	else:
            auxbas.append(map(int,line.split()[0:2])+map(float,line.split()[2:4]))
	    auxbas[i].insert(0,atom)

    return auxbas

# Calculate integrals of unnormalized auxiliary basis sets functions
# Corrected for S orbitals
# Pure functions for l>1 are assumed
def calc_aux_ints(auxbas):

    aux_ints = []
    for shell in auxbas:
        L = shell[2]
	a = shell[3]
	N = shell[4]
        if (L == 0): 
	    val = N*(math.pi/a)**(3./2.)
	else:
	    val = 0.0

	for i in range(2*L+1):
            aux_ints.append(val)

    return np.array(aux_ints)

# Calculate (P|r|Q) integrals
# Pure functions for l>1 are assumed
def calc_prq_ints(auxbas):

    prq_ints = []

    return np.array(prq_ints)
# read terms (P|Q)^(-1/2) where P,Q are indices of auxiliary (RI) basis functions
# and the integral contains the electronic repulsion
def read_2c2e_invsq_ints(fname, dims):

    f=open(fname)
    
    auxsqinv = np.zeros((dims,dims))
    for p in range(dims):
        for q in range(dims):
            line = f.readline().rstrip().split()
            auxsqinv[p,q] = float(line[0])

    f.close()

    return auxsqinv


# read terms (ij|P) where i,j are MO indices and P is auxiliary (RI) basis index
# and the integral involves the electronic repulsion
def read_3c2e_ints(fname, auxdims, no, nv):

    f=open(fname)

    auxmomo = np.zeros((no, nv, auxdims))
    for i in range(no):
        for p in range(auxdims):
            for j in range(nv):
                line = f.readline().rstrip().split()
                auxmomo[i,j,p] = float(line[0])

    f.close()

    return auxmomo


# read 2e-repulsion integrals in atomic orbitals; normally all are written out
# the Qchem output is only reliable for CARTESIAN basis sets
def read_ao_ints(fname, aodims):
 
    f=open(fname)
    ints = np.zeros((aodims, aodims, aodims, aodims))
    for i in range(aodims):
        for j in range(aodims):
            for k in range(aodims):
                for l in range(aodims):
                    line=f.readline().rstrip().split()
                    ints[i,j,k,l] = float(line[0])
    f.close()
    return ints


# read 2e-repulsion integrals in molecular orbitals
# only non-zero terms are written out in file 'fmosname'; 
# the corresponding indices are in different file 'findsname'
def read_mo_ints(fmosname, findsname, modims):

    finds=open(findsname)
    ints = np.zeros((modims, modims, modims, modims))
    line = finds.readline()

    with open(fmosname) as fmos:

        for line in fmos:
            inds = finds.readline()
            a,b,c,d = [(int(inds.rstrip().split()[i])-1) for i in range(4)] 
            ints[a,b,c,d] = float(line.rstrip().split()[0])
            ints[b,a,c,d] = ints[a,b,c,d]
            ints[a,b,d,c] = ints[a,b,c,d]
            ints[b,a,d,c] = ints[a,b,c,d]
            ints[c,d,a,b] = ints[a,b,c,d]
            ints[d,c,a,b] = ints[a,b,c,d]
            ints[c,d,b,a] = ints[a,b,c,d]
            ints[d,c,b,a] = ints[a,b,c,d]
    finds.close()
    return ints


# read 1e core integrals in atomic orbitals from Qchem output file
# core integrals are sum of kinetic energy and nuclear attraction integrals
def read_core_ints(fname, modims):

    core = np.zeros((modims, modims))
    f= open(fname)

    line = f.readline()
    while (not 'Core Hamiltonian Matrix' in line):
        line = f.readline()

    for i in range(0,modims/7*7,7):
        f.readline()
        for j in range(modims):
#            line = f.readline().rstrip().split()
            tmpline = f.readline()
            line=[]; line.append(tmpline[:5])
            line=line+re.findall(r"[-+]?\d{0,5}\.\d{5}",tmpline[5:])
            core[(int(line[0])-1),(i+0):(i+7)] = [float(line[k]) for k in range(1,8)]
            
    if(modims % 7):
        i = modims/7 * 7
        rest = modims % 7
        f.readline()
        for j in range(modims):
            line = f.readline().rstrip().split()
            core[(int(line[0])-1),(i+0):(i+rest)] = [float(line[k]) for k in range(1,rest+1)]
        
    f.close()
    return core


# read 1e nuclear integrals in atomic orbitals from Qchem output file
def read_nuc_ints(fname, modims):


    nuc = np.zeros((modims, modims))
    f= open(fname)

    line = f.readline()	
    while (not 'Nuclear Attraction Matrix' in line):
        line = f.readline()
    
    for i in range(0,modims/7*7,7):
        f.readline()
        for j in range(modims):
#            line = f.readline().rstrip().split()
            tmpline = f.readline()
            line=[]; line.append(tmpline[:5])
            line=line+re.findall(r"[-+]?\d{0,5}\.\d{5}",tmpline[5:])
            nuc[(int(line[0])-1),(i+0):(i+7)] = [float(line[k]) for k in range(1,8)]

    if(modims % 7):
        i = modims/7 * 7
        rest = modims % 7
        f.readline()
        for j in range(modims):
            line = f.readline().rstrip().split()
            nuc[(int(line[0])-1),(i+0):(i+rest)] = [float(line[k]) for k in range(1,rest+1)]

    f.close()
    return nuc


# read 1e kinetic energy integrals in atomic orbitals from Qchem output file
def read_kin_ints(fname, modims):

    
    kin = np.zeros((modims, modims))
    f= open(fname)

    line = f.readline() 
    while (not 'Kinetic Energy Matrix' in line):
        line = f.readline()
    
    for i in range(0,modims/7*7,7):
        f.readline()
        for j in range(modims):
#            line = f.readline().rstrip().split()
            tmpline = f.readline()
            line=[]; line.append(tmpline[:5])
            line=line+re.findall(r"[-+]?\d{0,5}\.\d{5}",tmpline[5:])
            kin[(int(line[0])-1),(i+0):(i+7)] = [float(line[k]) for k in range(1,8)]

    if(modims % 7):
        i = modims/7 * 7
        rest = modims % 7
        f.readline()
        for j in range(modims):
            line = f.readline().rstrip().split()
            kin[(int(line[0])-1),(i+0):(i+rest)] = [float(line[k]) for k in range(1,rest+1)]
    
    f.close()
    return kin


# read 1e overlap integrals in atomic orbitals from Qchem output file
def read_overlap_ints(fname, modims):


    smat = np.zeros((modims, modims))
    f= open(fname)

    line = f.readline()
    while (not 'Overlap Matrix' in line):
        line = f.readline()

    for i in range(0,modims/5*5,5):
        f.readline()
        for j in range(modims):
            line = f.readline().rstrip().split()
            smat[(int(line[0])-1),(i+0):(i+5)] = [float(line[k]) for k in range(1,6)]

    if(modims % 5):
        i = modims/5 * 5
        rest = modims % 5
        f.readline()
        for j in range(modims):
            line = f.readline().rstrip().split()
            smat[(int(line[0])-1),(i+0):(i+rest)] = [float(line[k]) for k in range(1,rest+1)]

    f.close()
    return smat

# read dipole/quadrapole matrix in atomic orbitals from Qchem output file
def read_dip_ints(fname, modims, x,y,z):

    dip = np.zeros((modims, modims))
    f= open(fname)

    line = f.readline()
    while (not 'Multipole Matrix (%s,%s,%s)' % (x,y,z) in line):
        line = f.readline()

    for i in range(0,modims/7*7,7):
        f.readline()
        for j in range(modims):
            line = f.readline().rstrip().split()
            dip[(int(line[0])-1),(i+0):(i+7)] = [float(line[k]) for k in range(1,8)]

    if(modims % 7):
        i = modims/7 * 7
        rest = modims % 7
        f.readline()
        for j in range(modims):
            line = f.readline().rstrip().split()
            dip[(int(line[0])-1),(i+0):(i+rest)] = [float(line[k]) for k in range(1,rest+1)]

    f.close()
    return dip


# read MO coefficients from file
def read_mo_coeff(fname, modims):

    f=open(fname)
    coeff = np.zeros((modims,modims))
    for i in range(modims):
        for j in range(modims):
            line = f.readline().rstrip().split()
            coeff[i,j] = float(line[0])
    f.close()
    return coeff


# transform a set of 2e integrals with a given set of transformation coefficients
def transform_2e(ints, coeff, dims):

    tints = np.zeros((dims, dims, dims, dims))
    for i in range(dims):
        for a in range(dims):
            tints[i,:,:,:] += ints[a,:,:,:] * coeff[i,a]

    tmp = np.zeros((dims, dims, dims, dims))
    for j in range(dims):
        for b in range(dims):
            tmp[:,j,:,:] += tints[:,b,:,:] * coeff[j,b]

    tints = np.zeros((dims, dims, dims, dims))
    for k in range(dims):
        for c in range(dims):
            tints[:,:,k,:] += tmp[:,:,c,:] * coeff[k,c]

    tmp = np.zeros((dims, dims, dims, dims))
    for l in range(dims):
        for d in range(dims):
            tmp[:,:,:,l] += tints[:,:,:,d] * coeff[l,d]

    return tmp


# transform a set of 1e integrals with a given set of transformation coefficients
def transform_1e(ints, coeff, dims):

    tints = np.zeros((dims, dims))
    for i in range(dims):
        for a in range(dims):
            tints[i,:] += ints[a,:] * coeff[i,a]

    tmp = np.zeros((dims, dims))
    for l in range(dims):
        for d in range(dims):
            tmp[:,l] += tints[:,d] * coeff[l,d]

    return tmp


# read a 2-index matrix to file
def read_2d_matrix(fname, M, N):

    f=open(fname,'r')
    mat=np.zeros((M,N))

    lines=f.read().splitlines()
    for i in range(len(lines)):
        tmp=lines[i].split()
        mat[int(tmp[0])-1,int(tmp[1])-1]=float(tmp[2])

    return mat


# print a 2-index matrix to file; only print terms bigger than given value
def print_2d_matrix(fname, mat, cutoff=0.):

    dims = mat.shape
    f = open(fname, 'w')

    for i in range(dims[0]):
        for j in range(dims[1]):
            if(abs(mat[i,j]) >= cutoff):
                #f.write(str(mat[i,j]))
                f.write("% .16e" % mat[i,j])
                f.write('\n')
    f.close()
    return

# print a 3-index matrix to file; only print terms bigger than given value
def print_3d_matrix(fname, mat, cutoff=1.0e-9):

    dims = mat.shape
    f = open(fname,'w')

    for i in range(dims[0]):
        for j in range(dims[1]):
            for k in range(dims[2]):
                if(abs(mat[i,j,k]) > cutoff):
                    f.write(str(i+1)+' '+str(j+1)+' '+str(k+1)+' '+str(mat[i,j,k]))
                    f.write('\n')
    f.close()
    return


# print a 4-index matrix to file; only print terms bigger than given value
def print_4d_matrix(fname, mat, cutoff=1.0e-9):

    dims = mat.shape
    f = open(fname,'w')

    for i in range(dims[0]):
        for j in range(dims[1]):
            for k in range(dims[2]):
                for l in range(dims[3]):
                    if(abs(mat[i,j,k,l]) > cutoff):
                        f.write(str(i+1)+' '+str(j+1)+' '+str(k+1)+' '+str(l+1)+' '+str(mat[i,j,k,l]))
                        f.write('\n')
    f.close()
    return

# print (ij|P) integrals
def print_ijP(auxmomo, fname='363.1', cutoff=1.0e-9):

    f= open(fname, 'w')

    for i in range(auxmomo.shape[0]):
        for j in range(auxmomo.shape[1]):
            for p in range(auxmomo.shape[2]):
                if(abs(auxmomo[i,j,p]) > cutoff):
                    f.write(str(i+1)+' '+str(j+1)+' '+str(p+1)+' '+str(auxmomo[i,j,p]))
                    f.write('\n')

    f.close()

    return 

# Get (P,Q)^(-1) 
def matrix_square(invsq):

    PQ_inv = np.dot(invsq,invsq)
    return PQ_inv


# multiply over auxiliary basis          
# input is indexed as(i,j,P) and (P,Q)
# output is indexed as (i,j,P)
def mult_ijP_PQ(auxmomo, invsq):

    out = np.tensordot(auxmomo, invsq, (2,0))
    return out

# multiply over auxiliary basis
# input is indexed as(i,j,P) and (i,j,P)
# output is indexed as (i,j,P)
# multiply matrices (ij|P)*(P|Q)^(-1/2) together to obtain approximated integrals (ij | kl)
# (ij | kl) ;= (ij|P)*(P|Q)^(-1/2) * (Q|R)^(-1/2)*(R|kl)
def mult_ijP_Pkl(auxmomo):

    two = np.tensordot(auxmomo, auxmomo, (2,2))
    return two

# make terms (ij|P)*(P|Q)^(-1/2) from scratch
def ri_sqrt_2eints(modims, auxdims, f_ijP = '362.0.txt', f_PQ = '363.0.txt'):

    invsq = read_2c2e_invsq_ints(f_PQ, auxdims)
    auxmomo = read_3c2e_ints(f_ijP, auxdims, modims, modims)
    vri = mult_ijP_PQ(auxmomo, invsq)
    return vri

# make approximate RI 2-e integrals from scratch
def ri_2eints(modims, auxdims, f_ijP = '362.0.txt', f_PQ = '363.0.txt'):

    vri = ri_sqrt_2eints(modims, auxdims, f_ijP, f_PQ)
    two = mult_ijP_Pkl(vri)    
    return two



