# Copyright (C) Isaure Chauvot de Beauchene (CNRS))

import sys, os
import numpy as np
from rmsdlib import multifit

def npy2to3(npy):
    if len(npy.shape) == 2:
        if npy.shape[1] == 3:
            npy = npy.reshape(1, npy.shape[0], npy.shape[1])
        else:
            npy = npy.reshape(npy.shape[0], int(npy.shape[1]/3), 3)
    else:
        assert len(npy.shape) == 3
    return npy

def npy3to2(npy):
    if len(npy.shape) == 3:
        npy = npy.reshape(npy.shape[0], 3*npy.shape[1])
    else:
        assert len(npy.shape) == 2 and npy.shape[1]%3 == 0
    return npy

def fit_multi_npy(a, ref):
    a = npy2to3(a)
    rotation, translation, RMSD = multifit(a, ref)
    rot = np.transpose(rotation, axes=(0,2,1))
    COM = a.sum(axis=1)/a.shape[1]
    centered = a - COM[:,None,:]
    rotated = np.einsum('...ij,...jk->...ik',centered,rot)
    fitted = rotated + COM[:,None,:]
    translated = fitted - translation[:,None,:]
    return translated, RMSD

def rmsdnpy(chains, pdb):
    chains = npy3to2(chains)
    reference = [ l for l in open(pdb).readlines() if l.startswith("ATOM") ]
    r = [ [float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in reference]
    ref = np.array(r)
    ref = ref.reshape(ref.shape[0]*ref.shape[1])
    ncoord = np.shape(chasins[0])[0]
    RMSD = [ (sum([(chain[i]-ref[i])**2 for i in range(ncoord)]) /(ncoord/3))**0.5 for chain in chains ]
    return RMSD

def map_npz(npz_file):
    print("map_npz",file=sys.stderr)
    sys.stderr.flush()
    npz = np.load(npz_file)
    nfrags = npz["nfrags"]
    poses, interactions =  [], []
    for n in range(nfrags-1):
        inter = npz["interactions-%d"%n]
        interactions.append(inter)
        poses.append(np.unique(inter[:,0]))
    poses.append(np.unique(inter[:,1]))
    npz = []
    #interactions = [ np.array(i, dtype=int) for i in inter]
    return interactions, poses
