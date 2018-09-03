# Given an xyz file, this script converts it into the q-chem
# input file format


import sys


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: xyz")
        sys.exit(1)

    xyz = sys.argv[1]
    prefix = xyz.split("/")[-1].split(".")[0]

    out = "$molecule\n"
    if "-" in prefix:
        out += "-1 1\n"
    else:
        out += "0 1\n"
    with open(xyz, "r") as f:
        f.readline()
        f.readline()
        for line in f:
            out += line
    out += "$end\n"

    sys.stdout.write(out)
