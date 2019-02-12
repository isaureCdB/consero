#!/usr/bin/env python3
# Copyright (C) Isaure Chauvot de Beauchene (CNRS))

import sys, os, argparse
import numpy as np
usage = "USAGE: python deredundant-npy.py structures.npy cutoff [--fit] [--center]\n"
from npy import npy3to2

############
parser =argparse.ArgumentParser(description=__doc__,
                        formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('npyfile', help="np array")
parser.add_argument('cutoff', help="cutoff",type=float)
parser.add_argument('--outp', help="np array deredundanted +/- fitted")
parser.add_argument("--fit",help="fit on 1st structure", action="store_true")
parser.add_argument("--center",help="center at origine", action="store_true")
parser.add_argument("--chunksize", type=int, default=1000)
args = parser.parse_args()
############

npy = np.load(args.npyfile)

name = args.npyfile.split(".npy")[0]
if not args.outp:
    outp = name + '-dr' + str(args.cutoff) + '.npy'
else:
    outp = args.outp

if args.fit:
    print("fitting", file=sys.stderr)
    from rmsdlib import multifit
    rotation, translation, RMSD = multifit(npy, npy[0])
    rot = np.transpose(rotation, axes=(0,2,1))
    COM = npy.sum(axis=1)/npy.shape[1]
    translated = npy - COM[:,None,:]
    #print(rot.shape, file=sys.stderr)
    #print(translated.shape, file=sys.stderr)
    npy = np.einsum('...ij,...jk->...ik',translated,rot)
    npy = npy[ [i for i in range(len(RMSD)) if RMSD[i] > args.cutoff or i==0] ]
    if not args.center:
        COM1 = npy[0].sum(axis=0)
        npy = npy + COM1

from cluster import cluster
nn = npy3to2(npy)
c = cluster(npy3to2(npy), args.cutoff, args.chunksize)
centers = np.array([ c[k][0] for k in list(c.keys())])
npclust = npy[centers]
clustlist = open('%s-clust%.1f.list'%(name,args.cutoff), "w")
print(name+"-clust"+str(args.cutoff))
for k in centers:
    print("%i "%(k+1), file=clustlist)
np.save(outp,npclust)
clustlist.close()
