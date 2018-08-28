# mpe-mol
Many-pair expansion for molecules

### Pre-requesite

Make sure you have libxc installed, which can be found [here](http://www.tddft.org/programs/libxc/).

### Compile

In the `src` directory, open `Makefile` with your favorite text editor, and set `LIBXC_ROOT` to the root directory of your libxc library. Then type
```
make -j8
```
which will generate a binary executable called `mpe.out`.

### Run

`mpe.out` can be launched in any directory by
```
${MPE_ROOT}/mpe.out alpha data_path
```
where `MPE_ROOT` is the path of the MPE `src` directory, `alpha` controls the strength of density localization, and `data_path` is the (relative or absolute) path where input data are stored.

### Preparing input data

Currently MPE works with data generated from Q-Chem. Several python+shell scripts that help generate these data can be found in the `script` directory. You can also find example input files in the `example` directory. Here we go through one of them for clarity.

* Step 1: get integrals from Q-Chem

   In `script`, type
   ```
   bash runqchem.sh ../example/bh/lda bh_lda
   ```
   which will 
   1. run the two Q-Chem input files, `bh_lda.in` and `bh_lda_grid.in`, in the `example/bh/lda` directory, and
   2. run two python scripts, `writedata.py` and `writegrid.py`, to grep appropriate data from Q-Chem output files.
   
* Step 2: determine xc functionals

   In addition to the integral files, you also need to tell MPE which exchange-correlation (xc) functional to use, which is controlled by file `xc_func.txt` in the same directory. It has two number in same row:
   ```
   func Kxx
   ```
   where `func` is an integer indexing different functionals, and `Kxx` is a portion of exact exchange (EXX) for global hybrid functional (of course, `0 < Kxx < 1`). Currently we only support a small number of functionals, and they are
   
   | `func` | functional name |
   |--------|-----------------|
   | 1      |    LDA          |
   | 2      |    PBE          |
   | 3      |    PBE0         |
   | 4      |    BLYP         |
   | 13     |    PBE0*        |
   
   For `func < 10`, the `Kxx` value in the input file will be overwritten by the appropriate value in the corresponding functional. For `func > 10`, you are free to tune `Kxx` for your own purpose.
