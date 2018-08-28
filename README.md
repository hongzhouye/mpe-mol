# mpe-mol
Many-pair expansion for molecules

### Pre-requesite

Make sure you have libxc installed, which can be found [here](http://www.tddft.org/programs/libxc/).

### Compile

In the `src` directory, open `Makefile` with your favorite editor, and set `LIBXC_ROOT` to the root directory of your libxc library. Then type
```
make -j8
```
which will generate a binary executable called `mpe.out`.

### Run

`mpe.out` can be launched in any directory by
```
${MPE_ROOT}/mpe.out alpha data_path
```
where `MPE_ROOT` is the path to the MPE `src` directory, `alpha` controls the strength of density localization, and `data_path` is the (relative or absolute) path where input data are stored.
