#!/usr/bin/env python3
# Copyright (C) Isaure Chauvot de Beauchene (CNRS)

import sys, os, argparse
import numpy as np

############
parser =argparse.ArgumentParser(description=__doc__,
                        formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--chainfile', help="ex: frag1-2_2A-meanrk1000.chains")
parser.add_argument('--outp', help="ex: frag1-2_2A-top1000.npy")
#parser.add_argument('--top', help="number of top-ranked chain to consider", action="store_true")
parser.add_argument('--npyfiles', nargs='+', help="frag1.npy frag2.npy...")
parser.add_argument('--preatoms', nargs='+', help="frag1.preatoms frag2.preatoms...")
parser.add_argument('--postatoms', nargs='+', help="frag1.postatoms frag2.postatoms...")
args = parser.parse_args()
############

assert len(args.preatoms) == len(args.postatoms)
assert len(args.preatoms) == len(args.npyfiles)

ll = [ l for l in open(args.chainfile).readlines()][1:]
if len(ll) < 2:  sys.exit()
nfrag = int( 0.5*(len(ll[0].split())- 4))

chains = [ [ int(j)-1 for j in l.split()[4:4+nfrag] ] for l in ll ]
Nchains = len(chains)

npy = [np.load(f) for f in args.npyfiles]
preatoms = [ [ int(l.strip())-1 for l in open(p).readlines() ] for p in args.preatoms ]
postatoms = [ [ int(l.strip())-1 for l in open(p).readlines() ] for p in args.postatoms ]

preat = []
for f in range(nfrag):
    indices = [ j for j in preatoms[f] if j not in postatoms[f] ]
    preat.append(npy[f][:,indices]) # npy[frag][chains][atoms]

midat = []
for f in range(nfrag):
    indices = [ j for j in preatoms[f] if j in postatoms[f] ]
    midat.append(npy[f][:,indices])

postat = []
for f in range(nfrag):
    indices = [ j for j in postatoms[f] if j not in preatoms[f] ]
    postat.append(npy[f][:,indices])

pdbchains = []

# c = pose of frag1, pose of frag2, pose of frag3
for c in chains:
    pdbchain = []
    pdbchain.extend(postat[0][c[0]])

    mix = 0.5*(midat[0][c[0]] + postat[1][c[1]])
    pdbchain.extend(mix)
    #print "b", len(pdbchain)
    #print np.array(pdbchain).shape, mix.shape

    for frag in range(nfrag-2):
        mix = (preat[frag][c[frag]] + midat[frag+1][c[frag+1]] + postat[frag+2][c[frag+2]]) /3
        pdbchain.extend(mix)
    #print "c", len(pdbchain)

    mix = 0.5*(preat[-2][c[-2]] + midat[-1][c[-1]])
    pdbchain.extend(mix)
    #print "d", len(pdbchain)

    pdbchain.extend(preat[-1][c[-1]])
    #print "e", len(pdbchain)

    pdbchains.append(pdbchain)

##############################
pdbchains = np.array(pdbchains)
#print pdbchains.shape
np.save(args.outp, pdbchains)

#print >> sys.stderr, np.shape(pdbchains)
#for i in range(np.shape(pdbchains[0])[0]/3):
#    #print >> sys.stdout, pdbchains[0][3*i], pdbchains[0][3*i+1], pdbchains[0][3*i+2]
