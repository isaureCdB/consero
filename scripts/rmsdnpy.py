#!/usr/bin/env python3
## Copyright (C) Isaure Chauvot de Beauchene (CNRS)


import sys, os, argparse
import numpy as np
from npy import npy2to3, npy3to2

def rmsdnpy(npy, refs, selatoms):
    struc = npy[:,selatoms,:]
    refs = refs[:,selatoms,:]
    nat = len(selatoms)
    a = struc[:,None,:,:]
    b = refs[None,:,:,:]
    d = a - b
    SD = np.einsum("ijkl, ijkl -> ij", d, d)
    print((np.shape(struc), np.shape(refs)), nat, file=sys.stderr)
    RMSD = np.sqrt(SD/nat)
    return RMSD

#######################################
parser =argparse.ArgumentParser(description=__doc__,
                        formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('npy', help="numpy coordinates")
parser.add_argument("--ref", help="reference pdb structures", nargs='+')
parser.add_argument("--sym", default=None, help="indices of symetrical atoms")
parser.add_argument("--atoms", default=None, help="indices of atoms to consider")
args = parser.parse_args()
#######################################

refs = args.ref

# convert references into lists of coordinates
if refs[0].split(".")[-1] == "npy":
    refs = [npy3to2(np.load(r)) for r in refs]
else:
    coors = []
    for r in refs:
        ll = [ l for l in open(r).readlines() if l.startswith("ATOM") ]
        coor = [ [float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in ll]
        coors.append(coor)
    refs = np.array(coors)

npy = npy2to3(np.load(args.npy))
nat = np.shape(npy)[1]
selatoms = range(nat)
if args.atoms:
    selatoms = [int(i)-1 for i in args.atoms]

rmsds = rmsdnpy(npy, refs, selatoms)
rmsd = np.min(rmsds, axis=1)
print(rmsd.shape)

for i in range(len(rmsd)):
    print('%i %.2f'%(i+1, rmsd[i]))
