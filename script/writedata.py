from readints import *
from sys import argv
#from HF import *
import numpy as np
import time
import os

def writedata(name, rd, cwd, fcore="802.0.txt", fints='22.0.txt', finds = '809.0.txt', fmos='53.0.txt', f_PQ='363.0.txt', f_ijP='362.0.txt', cutoff=0):

    # read number of basis functions in Qchem output file
    # modims, auxdims = read_number_mos(name+'.out')
    #@@HY: the function above is too slow for giant output file!
    #      I wrote an efficient implementation which 
    #      (i) greps nbas/nelec/naux much faster and
    #      (ii) extracts Ekin and Enuc matrices so no need to read
    #           them from the output file line by line (extremly slow!)
    os.system("bash {:s}/grep_int_dim.sh {:s}/{:s}.out".format(rd, cwd, name))
    with open("{:s}/dimensions.txt".format(cwd), "r") as f:
       modims, nelec, auxdims = list( map(int, f.readline().split()))

    # atom coordinates
#    natom = read_atom_number(name+'.out')
#    coor = read_atom_coordinate(name+'.out',natom)

    # MO coefficients^M
    print("Reading MO coeff...", flush=True)
    coeff = read_mo_coeff(fmos, modims)
#    print_2d_matrix(name+'_mos.txt', coeff, cutoff)

    # 1-electron (core = kinetic energy + nuclear attraction) integrals
    # transfrom to MO basis
    # core = read_core_ints(name+'.out', modims)
    # @@HY: more accurate core Hamiltonian from binary files
    print("Reading core Hamiltonian...", flush=True)
    core = np.loadtxt(fcore).reshape(modims, modims)
    oneeint = transform_1e(core, coeff, modims)

    print("Reading core Hamiltonian component: kinetic energy...", flush=True)
    # kin = read_kin_ints(name+'.out', modims)
    kin = read_kin_ints("_hy_tmp_kin", modims)
    oneekin = transform_1e(kin, coeff, modims)     

    print("Reading core Hamiltonian component: nuclear attraction...", flush=True)
    # nuc = read_nuc_ints(name+'.out', modims)
    nuc = read_nuc_ints("_hy_tmp_nuc", modims)
    oneenuc = transform_1e(nuc, coeff, modims)

    # @@HY: haven't seen this been used in anywhere so commented
    # print("Reading AO overlap matrix...", flush=True)
    # overlap = read_overlap_ints(name+'.out', modims)

    # coefficients (P | Q)^(-1/2)
    print("Reading aux basis overlap matrix (P|Q)...", flush=True)
    invsq = read_2c2e_invsq_ints(f_PQ, auxdims)
#    print_2d_matrix(name+'_invsqPQints.txt', invsq, cutoff)

    # coefficients (ij | P)
    print("Reading aux basis-MO orbital overlap matrix (ij|P)...", flush=True)
    auxmomo = read_3c2e_ints(f_ijP, auxdims, modims, modims)
#    print_3d_matrix(name+'_ijPints.txt', auxmomo, cutoff)

    # @@HY: just realize the 2e integral read in was never used...commented
    # 2e integrals
    # print("Reading 2e Hamiltonian...", flush=True)
    # twoeint = read_mo_ints(fints, finds, modims)
#    print_4d_matrix(name+'_2eints.txt', twoeint, cutoff)
#    print '2e integral',twoeint

    # calculate (P | Q)^(-1)
    print("Computing (P|Q)...", flush=True)
    PQ_inv = matrix_square(invsq)

    # calculate (P | Q)    
    PQ = np.linalg.inv(PQ_inv)

    # read the number of occupied orbitals
    nocc = read_number_nocc(name+'.out')

    # read auxiliary(RI) basis set    
#    auxbas = read_aux_basis(name+'.out')
#    print auxbas

    # Represent the density in RI basis
#    aux_coeff = dens_aux_coeff(nocc,auxdims,PQ_inv,auxmomo)
#    print 'aux_coeff_alpha',aux_coeff/2.

    # Integrate each RI basis function
#    aux_ints = calc_aux_ints(auxbas)

    # Calculate the total density in RI basis 
#    dens_aux = np.dot(aux_coeff,aux_ints)
#    print 'dens_aux',dens_aux

    # read multipole matrix and transfrom to MO basis
#    dipx = read_dip_ints(name+'.out', modims,1,0,0)
#    dipx_mo = transform_1e(dipx, coeff, modims)

#    dipy = read_dip_ints(name+'.out', modims,0,1,0)
#    dipy_mo = transform_1e(dipy, coeff, modims)

#    dipz = read_dip_ints(name+'.out', modims,0,0,1)
#    dipz_mo = transform_1e(dipz, coeff, modims)

#    dipx2 = read_dip_ints(name+'.out', modims,2,0,0)
#    dipx2_mo = transform_1e(dipx2, coeff, modims)

#    dipy2 = read_dip_ints(name+'.out', modims,0,2,0)
#    dipy2_mo = transform_1e(dipy2, coeff, modims)

#    dipz2 = read_dip_ints(name+'.out', modims,0,0,2)
#    dipz2_mo = transform_1e(dipz2, coeff, modims)

    print("Dumping...", flush=True)
    # @@HY: "dimensions.txt" has been taken care of by "grep_int_dim.sh"
    # f1=open('dimensions.txt','w')
    # f1.write(str(modims)+' '+str(nocc)+' '+str(auxdims))
    # f1.close()

    print_2d_matrix('oneeint.txt', oneeint, cutoff)
    print_2d_matrix('oneekin.txt', oneekin, cutoff)
    print_2d_matrix('oneenuc.txt', oneenuc, cutoff)
    print_2d_matrix('PQ.txt', PQ, cutoff)
    print_2d_matrix('PQ_inv.txt', PQ_inv, cutoff)
    print_2d_matrix('coeff.txt', coeff, cutoff)
    # print_2d_matrix('overlap.txt',overlap,cutoff)

start_time = time.time()
name = argv[1]

# root directory where this script is stored
rd = os.path.dirname(os.path.realpath(__file__))
# current working directory
cwd = os.getcwd()

writedata(name, rd, cwd)
os.system("rm {:s}/_hy_tmp_kin {:s}/_hy_tmp_nuc".format(cwd, cwd))
print("--- %s seconds ---" % (time.time() - start_time))
