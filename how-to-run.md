# 10 steps to run an MPE calculation

Take BH3 molecule as an example. Suppose everythin is run in some directory given by `$run_path`.

### Generate a Q-Chem input file for grid

```shell
cd $run_path
mkdir grid; cd grid
cat ${MPE_ROOT}/molecule/existing_cases/bh3.in >> bh3_grid.in    # molecule section
echo "" >> bh3_grid.in
cat ${MPE_ROOT}/rem/grid.in >> bh3_grid.in                       # rem section
echo "" >> bh3_grid.in
echo "\$basis" >> bh3_grid.in                                    # basis section
cat ${MPE_ROOT}/basis/uTZ_ordered/B.in >> bh3_grid.in            # basis for B
echo "****" >> bh3_grid.in
cat ${MPE_ROOT}/basis/uTZ_ordered/H.in >> bh3_grid.in            # basis for H
echo "****" >> bh3_grid.in
echo "\$end" >> bh3_grid.in
```

### Run it and get the grid

```shell
cd ${MPE_ROOT}/script
bash get_qchem_grid.sh $run_path/grid bh3                        # you might need to change the scratch directory in this script
```

You should see three files generated in `$run_path/grid`: `grid_dimension.txt`, `grid.txt`, `ao_grid.txt`.

### Generate another input file for MO/aux integrals

This step is functional-dependent. Take pbe as an example.

```shell
cd $run_path
mkdir data; cd data
cat ${MPE_ROOT}/molecule/existing_cases/bh3.in >> bh3_pbe.in    # molecule section
echo "" >> bh3_pbe.in
cat ${MPE_ROOT}/rem/scf.in >> bh3_pbe.in                        # rem section
echo "" >> bh3_pbe.in
cat ${MPE_ROOT}/xc/pbe.in >> bh3_pbe.in                         # xc_functional section
echo "" >> bh3_pbe.in
echo "\$basis" >> bh3_pbe.in                                    # basis section
cat ${MPE_ROOT}/basis/uTZ_ordered/B.in >> bh3_pbe.in            # basis for B
echo "****" >> bh3_pbe.in
cat ${MPE_ROOT}/basis/uTZ_ordered/H.in >> bh3_pbe.in            # basis for H
echo "****" >> bh3_pbe.in
echo "\$end" >> bh3_pbe.in
echo "" >> bh3_pbe.in
cat >> bh3_pbe.in << EOF                                        # second input file; for aux basis
@@@

\$molecule
    read
\$end

EOF
cat ${MPE_ROOT}/rem/aux.in >> bh3_pbe.in                        # rem section
echo "" >> bh3_pbe.in
echo "\$basis" >> bh3_pbe.in                                    # basis section
cat ${MPE_ROOT}/basis/uTZ_ordered/B.in >> bh3_pbe.in            # basis for B
echo "****" >> bh3_pbe.in
cat ${MPE_ROOT}/basis/uTZ_ordered/H.in >> bh3_pbe.in            # basis for H
echo "****" >> bh3_pbe.in
echo "\$end" >> bh3_pbe.in
echo "" >> bh3_pbe.in
echo "\$aux_basis" >> bh3_pbe.in                                # aux_basis section
cat ${MPE_ROOT}/basis/aux_uTZ/B.in >> bh3_pbe.in                # aux basis for B
echo "****" >> bh3_pbe.in
cat ${MPE_ROOT}/basis/aux_uTZ/H.in >> bh3_pbe.in                # aux basis for H
echo "****" >> bh3_pbe.in
echo "\$end" >> bh3_pbe.in
```

### Run it and get the needed data (i.e. integrals)

```shell
cd $MPE_ROOT/script
bash get_qchem_int.sh $run_path/data/pbe bh3_pbe                # you might need to change the scratch directory in this script
```
You should get a bunch of `txt` files including `oneeint.txt`, `oneekin.txt`, and etc.

### Doing SAH decomposition

```shell
cd $run_path
mkdir alp_4.0; cd alp_4.0                                       # use alpha=4.0 for the moment
mkdir pbe; cd pbe
$MPE_ROOT/src/mpe.out 4.0 ../../data/pbe ../../grid
```
This will crash after running the SAH decomposition with an error message like
```
At line 192 of file RunMPE.F (unit = 21, file = 'xc_func.txt')
Fortran runtime error: End of file
```
which is intended! The point is that it writes the pair densities information in `pairs.txt` which we can resuse in the next step.

### 
