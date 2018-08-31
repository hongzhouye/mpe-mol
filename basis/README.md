Collected here are basis sets for needed atoms, including

H, B, Be, O, N, F, Cl 

For each atom, the normal basis set to expand MOs are stored
in directory "uTZ"; the corresponding auxiliary basis which
is obtained by removing some of the most contracted and diffuse
functions from the original uTZ basis set is stored in "aux_uTZ".

To use them, collect bases for all atoms you want in one file,
with different atoms separated by a line containing four asterisk,
e.g.

H   0
S   1   1.00
...
****
O   0
S   1   1.00
...
****
N   0
S   1   1.00
...
