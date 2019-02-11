#!/usr/bin/env python3

import sys, os
import numpy as np
sys.path.append(os.environ["ATTRACTTOOLS"])
from rmsdlib import multifit
from npy import npy2to3, npy3to2

def rmsdnpy(npy, ref, selatoms):
    #print(npy.shape, file=sys.stderr)
    ref = ref.reshape(ref.shape[0]*ref.shape[1])
    print((np.shape(npy), np.shape(ref)), file=sys.stderr)
    struc = npy3to2(npy[:,selatoms,:])
    nat = len(selatoms)
    RMSD = [ (sum([(s[i]-ref[i])**2 for i in range(nat*3)]) /nat)**0.5 for s in struc ]
    return RMSD

def pairwise_rmsdnpy(a, b):
    a = npy3to2(a)
    b = npy3to2(b)
    ncoord = np.shape(a[0])[0]
    print((np.shape(a), np.shape(b)), file=sys.stderr)
    RMSD = np.array([[(sum([(a1[i]-b1[i])**2 for i in range(ncoord)]) /(ncoord/3))**0.5  for b1 in b] for a1 in a])
    return RMSD

#######################################
parser =argparse.ArgumentParser(description=__doc__,
                        formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('npy', help="numpy coordinates")
parser.add_argument("--ref", default=None, help="reference pdb structures", nargs='+')
parser.add_argument("--sym", default=None, help="indices of symetrical atoms")
args = parser.parse_args()
#######################################

npyfile = sys.argv[1]    # npy.npy
if sys.argv[2].split(".")[-1] == "npy":
    ref = np.load(sys.argv[2])
    if len(ref.shape) > 2:
        assert ref.shape[2] == 3
        ref = ref.reshape((ref.shape[1],3))
else:
    pdb = sys.argv[2] # Lboundr.pdb
    ll = [ l for l in open(pdb).readlines() if l.startswith("ATOM") ]
    r = [ [float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in ll]
    ref = np.array(r)

npy = npy2to3(np.load(npyfile))
nat = np.shape(npy)[1]
selatoms = range(nat)
if len(sys.argv) > 3:
    selatoms = [int(i)-1 for i in sys.argv[3:]]

RMSD = rmsdnpy(npy, ref, selatoms)
for i in range(len(RMSD)):
    print('%i %.2f'%(i+1, RMSD[i]))
