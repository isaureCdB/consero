## Copyright (C) Isaure Chauvot de Beauchene (CNRS)

#!/usr/bin/env python3
import numpy as np
import sys, argparse

########################
parser =argparse.ArgumentParser(description=__doc__,
formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('npyfile')
parser.add_argument('outpfile')
parser.add_argument('--atom', '--atoms', type=int, nargs='+', help="indices of atoms to select")
parser.add_argument('--structure', '--structures', '--struct', dest='structure',
                    type=int, nargs='+', help="indices of structures to select")
parser.add_argument('--atname', '--atnames','--atomname', dest='atname', type=str,
                    nargs='+', help="name of atoms to select")
parser.add_argument('--top', help="number of top structures to select", type=int)
parser.add_argument('--template', help="pdb template to select atom indices")

args = parser.parse_args()
########################

npy = np.load(args.npyfile)
reshape = False
if len(npy.shape) == 3:
    assert npy.shape[2] == 3
else:
    assert len(npy.shape) == 2 and npy.shape[1]%3 == 0
    npy = npy.reshape(npy.shape[0], int(npy.shape[1]/3),3)
    reshape = True

if args.structure:
    sel = [ int(i)-1 for i in args.structure]
    npy = npy[sel, : , :]

sel = range(npy.shape[1])
if args.atom:
    sel = [ int(i)-1 for i in args.atom]
elif args.atname:
    if args.template is None:
        print("atname option requires a pdb template (--template)", file=sys.stderr)
        raise
    pdb = [l for l in open(args.template).readlines() if l.startswith('ATOM') ]
    names = set(args.atname)
    sel = [ lnr for lnr, l in enumerate(pdb) if l[13:16].strip() in names ]

npy = npy[:, sel , :]

if reshape:
    npy = npy.reshape(npy.shape[0],npy.shape[1]*3)

np.save(args.outpfile, npy)

'''
npy = np.load(sys.argv[1])
if args.top is not None:
    sel=range(args.top)
else:
    try:
        sel = int(sys.argv[2]) - 1
        outp = "%s-%i.npy"%(sys.argv[1].split(".npy")[0], sel+1)
        if len(sys.argv) > 3:
            outp = sys.argv[3]
        print(sel, file=sys.stderr)
    except:
        selection = [ int(l.split()[0]) for l in open(sys.argv[2]).readlines() ]
        sel = [ i-1 for i in selection]
        outp = sys.argv[3]

reshape = False
if len(npy.shape) == 3:
    assert npy.shape[2] == 3
else:
    assert len(npy.shape) == 2 and npy.shape[1]%3 == 0
    npy = npy.reshape(npy.shape[0], int(npy.shape[1]/3),3)
    reshape = True

npy2 = npy[sel, : , :]

#if reshape:
#    npy2 = npy2.reshape(npy2.shape[0],npy2.shape[1]*3)

np.save(outp, npy2)
'''
