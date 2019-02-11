# Copyright (C) Isaure Chauvot de Beauchene & Sjoerd J. de Vries (TUM)

#!/usr/bin/env python2

import sys, json, argparse
import numpy as np
from math import log
from npy import npy2to3

###########
parser =argparse.ArgumentParser(description=__doc__,
                        formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--graph', help="number of fragments to assemble")
parser.add_argument('--meanrank', help="maximal geometric mean of pose ranks in a chain")
parser.add_argument('--preatoms', nargs='+', help="frag1-preatoms.npy frag2-preatoms.npy ...")
parser.add_argument('--postatoms', nargs='+', help="frag1-postatoms.npy frag2-postatoms.npy ...")
parser.add_argument('--rmsd', nargs='+', help="frag1.rmsd frag2.rmsd ...", default=None)

args = parser.parse_args()
############

tree = json.load(open(args.graph))
max_meanrank = float(args.meanrank)
nfrag = tree["nfrags"]

postatoms = [npy2to3(np.load(f)) for f in args.postatoms]
preatoms = [npy2to3(np.load(f)) for f in args.preatoms]

meanrank_threshold = log(max_meanrank) * nfrag

lrmsds = []
if args.rmsd is not None:
    for f in args.rmsd:
      lrmsd = {}
      lnr = 0
      for l in open(f):
        lnr += 1
        ll = l.split()
        if len(ll) != 2: continue
        k = ll[0]
        if k == "l-RMSD": k = lnr
        else: k = int(k)
        value = float(ll[1])
        lrmsd[k] = value
      lrmsds.append(lrmsd)


fragments = []
interactions = [{} for l in range(nfrag-1)]
for cnr, tclus in enumerate(tree["clusters"]):
  clus = []
  for tcc in tclus:
    assert tcc["radius"] == 0 #for now, limit ourselves to fully collapsed trees
    rank = np.array(tcc["ranks"], dtype="int")[0]
    clus.append(rank)
  fragments.append(clus)
for cnr, tinter in enumerate(tree["interactions"]):
  inter = interactions[cnr]
  for source, target in tinter:
    if source not in inter: inter[source] = []
    inter[source].append(target)

count = 0
def chain(indices, sum_overlap_msds, meanrank):
    global count
    lr = []
    if args.rmsd is not None:
        for inr, i in enumerate(indices):
            lr.append(lrmsds[inr][i])
    lrms = np.sqrt(sum([v*v for v in lr])/len(lr))
    o = np.sqrt(sum_overlap_msds / (len(indices) - 1) )
    meanrank = "%.3f" % meanrank
    lrms = "%.3f" % lrms
    o = "%.3f" % o
    print "#indices", lrms, meanrank, o,
    for i in indices: print i,
    for l in lr: print l,
    print
    count += 1

def walk(pos, curr, indices, sum_overlap_msds, curr_meanrank):
  ind = fragments[pos][curr]
  new_indices = indices + (ind,)
  curr_meanrank += log(ind)
  if curr_meanrank > meanrank_threshold: return
  overlap_msd = 0
  if pos > 0:
    pre = preatoms[pos-1][indices[-1]-1]
    post = postatoms[pos][ind-1]
    d = post - pre
    overlap_msd = (d*d).sum()/pre.shape[0]
  new_sum_overlap_msds = sum_overlap_msds + overlap_msd
  if pos == nfrag - 1:
    meanrank = np.exp(curr_meanrank/nfrag)
    chain(new_indices, new_sum_overlap_msds, meanrank)
    return
  for target in interactions[pos][curr]:
    walk(pos+1, target, new_indices, new_sum_overlap_msds, curr_meanrank)

print "#header <mean (root-mean-sq) ligand rmsd> <mean (geometric mean) rank>  <rms-overlap-rmsd> <indices> <ligand rmsds>"
f = fragments[0]
for ff in range(len(f)):
  walk(0, ff, (), 0, 0)

print >> sys.stderr, count
