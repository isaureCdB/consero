#!/usr/bin/env python3

import sys, argparse
import numpy as np
from npy import npy2to3
########################
parser =argparse.ArgumentParser(description=__doc__,
                        formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('npy')
parser.add_argument('template')
parser.add_argument('--list', help = "list of structure indexes")
parser.add_argument('--index', help = "structure indexe", type = int)
parser.add_argument('--insert', action = "store_true")
args = parser.parse_args()
########################
npy = npy2to3(np.load(args.npy))

def parse_template_insertions(template):
    template = []
    template_insertions = {}
    for l in open(args.template).readlines():
        if l.startswith("ATOM") or l.startswith("HETATM"):
            template.append(l)
        else:
            pos = len(template)
            if pos not in template_insertions:
                template_insertions[pos] = []
            template_insertions[pos].append(l.rstrip())
    return template, template_insertions

def parse_template(template):
    template = []
    for l in open(args.template).readlines():
        if l.startswith("ATOM") or l.startswith("HETATM"):
            template.append(l)
    return template

def convert_insertions(npy, template, template_insertions):
    assert len(template) == npy.shape[1], (len(template), npy.shape[1])
    m = 1
    Nat = len(template)
    for j in range(npy.shape[0]):
        print("MODEL "+str(m))
        for i in range(Nat):
            if i in template_insertions:
                for line in template_insertions[i]:
                    print(line)
            l = template[i]
            if len(l) < 54:
                print(l)
                continue
            coor = npy[j][i]
            print(l[:30]+"%8.3f%8.3f%8.3f"%(coor[0], coor[1], coor[2])+l[54:-1])
        m+=1
        print("ENDMDL")

def convert(npy, template):
    assert len(template) == npy.shape[1], (len(template), npy.shape)
    m = 1
    for j in range(npy.shape[0]):
        print("MODEL "+str(m))
        for nl, l in enumerate(template):
            coor = npy[j][nl]
            print(l[:30]+"%8.3f%8.3f%8.3f"%(coor[0], coor[1], coor[2])+l[54:-1])
        m+=1
        print("ENDMDL")

if args.list:
    sel = [ int(l.split()[0])-1 for l in open(args.list).readlines() ]
    npy = npy[sel,:,:]
if args.index:
    npy = npy[args.index-1:args.index,:,:]

if not args.insert:
    template = parse_template(args.template)
    convert(npy, template)
else:
    template, template_insertions = parse_template_insertions(args.template)
    convert_insertions(npy, template, template_insertions)
