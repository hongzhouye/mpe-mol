import numpy as np
import sys


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: raw_basis")
        sys.exit(1)

    inp = sys.argv[1]
    prefix = inp.split("/")[-1].split(".")[0]

    # determine # of lines (inefficient, but who cares)
    with open(inp, "r") as f:
        nline = 0
        for line in f:
            nline += 1

    out = "{:s}   0\n".format(prefix)
    with open(inp, "r") as f:
        read = 0
        iline = 0
        for line in f:
            if prefix in line:
                angmom = line.split()[1]
                read = 0
            else:
                out += "{:s}  1  1.00\n".format(angmom)
                exponent = float(line.split()[0])
                out += "{:17.7f}\t1.0000".format(exponent)
                read = 1

            iline += 1
            if read and iline != nline:
                out += "\n"

    print(out)
