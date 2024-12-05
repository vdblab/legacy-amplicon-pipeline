#!/usr/bin/env python

import sys
import os

def main(infile, outdir):
    print("reading input file")
    with open(infile, "r") as inf:
        for i, line in enumerate(inf):
            if i > 0:
                samp = line.split("\t")[0]
                print("  - " + samp)
                with open(os.path.join(outdir, samp + ".sample"), "w") as outf:
                    outf.write(samp + "\n")
    # add an unassigned one
    with open(os.path.join(outdir,"Unassigned.sample"), "w") as outf:
        outf.write("Unassigned" + "\n")


with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    print(snakemake.__dict__)
    print("\n")
    if not os.path.exists(snakemake.params[0]):
        print("creating dir: %s\n" % snakemake.params[0])
        os.mkdirs(snakemake.params[0])
    main(infile = snakemake.input[0], outdir = snakemake.params[0])
