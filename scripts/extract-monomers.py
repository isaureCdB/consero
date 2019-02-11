#!/usr/bin/env python3

import sys, os, argparse, numpy as np

############
parser =argparse.ArgumentParser(description=__doc__,
                        formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('pdbfiles', help="list of PDB files of docking poses")
args = parser.parse_args()
############

#motifs = [l.split() for l in open(args.resnames).readlines()]

counts = {}
coor = {}

def save_coor(reslines, newres):
    global counts, coor
    if not len(reslines):
        return
    if newres not in counts:
        counts[newres] = 0
        coor[newres] = []
    counts[newres] += 1
    c = []
    for l in reslines:
        c.append([ float(i) for i in [l[30:38], l[38:46], l[46:54]] ])
    coor[newres].append(c)

reslines = []
res = None
for pdb in open(args.pdbfiles):
    for l in open(pdb.strip()):
        if not l.startswith("ATOM"):
            res = None
            continue
        newres = l[17:20].strip()
        if newres != res:
            if res != None:
                save_coor(reslines, res)
            reslines = []
            res = newres
        reslines.append(l)

for res in coor:
    print(res)
    np.save("%s.npy"%res, np.array(coor[res]))

'''
ATOM      1  O5  SGN     1       1.465   1.228  51.133 5
'''
