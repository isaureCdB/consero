#!/usr/bin/env python2

import sys, os
sys.path.append(os.environ["ATTRACTTOOLS"])
from rmsdlib import multifit
import numpy as np
from npy import npy2to3, npy3to2

a = np.load(sys.argv[1])
cutoff = float(sys.argv[2])
output_npy = open(sys.argv[1].split(".npy")[0]+"-clust"+str(cutoff)+".npy","w")
output_list = open(sys.argv[1].split(".npy")[0]+"-clust"+str(cutoff),"w")

assert len(sys.argv)==3, "arg: file.npy cutoff"

a = npy2to3(a)
assert len(a.shape)==3 and a.shape[2] == 3
rmsd_matrix = np.zeros((len(a), len(a)))
for n in range(len(a)):
    rmsds = multifit(a, a[n])[2]
    rmsd_matrix[n] = rmsds

neighbor_matrix = (rmsd_matrix < cutoff)

clustered = 0
clusters = []
while clustered < len(neighbor_matrix):
    neigh = neighbor_matrix.sum(axis=0)
    heart = neigh.argmax()
    leaf = np.where(neighbor_matrix[heart])[0]
    for cs in leaf:
        neighbor_matrix[cs,:] = False
        neighbor_matrix[:, cs] = False
    new_clust = [heart+1] + [v+1 for v in leaf if v != heart]
    clusters.append(new_clust)
    clustered += len(new_clust)

for c in clusters:
    for i in c:
        print >> output_list, i,
    print >> output_list,""

clust = [c[0]-1 for c in clusters]
clust_struc = a[clust]
np.save(output_npy, clust_struc)
output_list.close()
