#!/usr/bin/env python3
# Copyright (C)  Isaure Chauvot de Beauchene (CNRS)

import numpy as np
import sys, argparse, os
from rmsdlib import multifit
from npy import npy2to3, fit_multi_npy

############
parser =argparse.ArgumentParser(description=__doc__,
                        formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('inp_npy')
parser.add_argument('outp_npy')
parser.add_argument('--lib', nargs='+', help="list of npy monomer libraries")
parser.add_argument('--noext', help="do not replace termini", action="store_true")
args = parser.parse_args()
############
chains = np.load(args.inp_npy)
libraries = [ np.load(lib) for lib in args.lib ]
Nres = len(libraries)
Nchains = len(chains)

chains = npy2to3(chains)

atoms = [0] # first atom in each residue
for lib in libraries:
    atoms.append(atoms[-1]+lib.shape[1])

monomers_chains = [ chains[:, atoms[i]:atoms[i+1] ,:] for i in range(Nres) ]

r = range(Nres)
if args.noext:
    r = r[1:-1]

for j in range(Nchains):
    for i in r:
        monomer = chains[j, atoms[i]:atoms[i+1],:]
        fitted, RMSD = fit_multi_npy(libraries[i], monomer)
        best = fitted[ RMSD.argmin(axis=0)]
        monomer[:] = best

print(chains.shape)
np.save(args.outp_npy, chains)
