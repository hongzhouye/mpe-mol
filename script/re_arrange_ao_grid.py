import numpy
import sys


def ao_grid_xform(f1, f2, fout):
    with open(f1, "r") as f:
        buf = list(map(int, f.readline().split()))
        natom = buf[0]
        assert(len(buf) == natom+2)
        nbas = buf[-1]
        np_per_atom = buf[1:-1]

    np = sum(np_per_atom)
    np0_at_atom = [sum(np_per_atom[:i]) for i in range(natom+1)]
    print(np0_at_atom)
        
    grid = numpy.zeros([4, np*nbas])
    print(grid.shape)
    with open(f2, "r") as f:
        for i in range(natom):
            np0 = np_per_atom[i]
            for j in range(nbas):
                block = numpy.array([list(map(float, f.readline().split())) \
                    for k in range(np0)]).T
                grid[:, np0_at_atom[i]+j*np:np0_at_atom[i+1]+j*np] = block
                # for debug
                # print(block)
                # print("---")
                # print(grid)
                # print("===================")

    numpy.savetxt(fout, grid, fmt="% .8f")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: grid_dimension.txt ao_grid.txt new_ao_grid.txt")
        sys.exit(1)

    ao_grid_xform(sys.argv[1], sys.argv[2], sys.argv[3])
