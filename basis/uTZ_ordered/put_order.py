import sys
import numpy as np


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: inp_basis")
        sys.exit(1)

    inp = sys.argv[1]
    
    out = ""
    buf = []
    with open(inp, "r") as f:
        out += f.readline()
        buf = [line for line in f]

    n = len(buf)

    # grep exponents
    exps = [float(b.split()[0]) for b in buf[1:n:2]]

    # order buf
    for count, j in enumerate(np.argsort(exps)[::-1]):
        out += buf[2*j]    
        if count == n//2-1:
            out += buf[2*j+1].replace("\n", "")
        else:
            out += buf[2*j+1]
    
    print(out)
