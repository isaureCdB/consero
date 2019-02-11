#!/usr/bin/env python3

import sys
from collections import defaultdict

def dict_from_template(template):
    dico = defaultdict(str)
    resname = None
    a = []
    for l in template:
        if resname not in list(dico.keys()) and l[17:20] == resname:
            a.append(l[13:16])
        if l[17:20] != resname or l == template[-1]:
            if len(a) > 1:
                dico[resname] = a
            a = [l[13:16]]
            resname = l[17:20]
    return dico

def invdict(inpdict):
    invdict = defaultdict(str)
    for a in list(inpdict.keys()):
        for b in inpdict[a]:
            invdict[b] = a
    return invdict

def printmodel(atindices, residues, nres, struc):
    j=0
    print("MODEL "+str(struc))
    for i in range(nres):
        for l in [l for (a,l) in sorted(zip(atindices[i], residues[i])) ]:
            j+=1
            print("ATOM" + (7-len(str(j)))*" " + str(j) + " " + l[11:54]+l[-1])
    print("ENDMDL")

# pdbfile = Docking poses in PDB format.
# can be a sinle- or multi-pdb file.
pdbfile = sys.argv[1]
pdb = [l for l in open(pdbfile).readlines()] #input

# indexfile = file with one line per residue, each line i containing the
# indices (from 1) of the atoms in the input pdb belonging to residue i.
#
# If the file contains more indices than atoms in the pdbfile, the additional
# indices are ignored
indexfile = sys.argv[2]
indices = [ [int(k) for k in l.split()] for l in open(indexfile).readlines()]

# Template PDB file for the order of the atom in the output files.
# Can be any coordinates. Only the atom's naming is read out.
# The order of the atoms in the output is determined by the 'indices' input file :
# Atom of indice x in the indices file will be at position x in the outpout file,
# and will be given the name of the nth atom in the template PDB.
#
templatefile = sys.argv[3]
template = [l for l in open(templatefile).readlines() if len(l) > 13 ]

nres, natoms = len(indices), sum( [len(i) for i in indices] )

dict_res2at = defaultdict(str)
for i in range(len(indices)):
    dict_res2at[i] = indices[i]

dict_inv_names = {
            "SGN": ["SGN", "47Y", "07Y", "49Y"],
            "IDS": ["IDS", "UAP", "42u", "02u", "42U", "02U"],
		    "BDP" : ["4ZB", "0ZB", "BDP"],
		    "ASG" : ["34V", "04V", "ASG"],
		    "NAG" : ["3YB", "NAG", "0YB"],
            "AHR" : ["AHR"],
			}

dict_names = invdict(dict_inv_names)
dict_atoms =  dict_from_template(template)
dict_at2res = invdict(dict_res2at)
dict_atname = { "O1S":"OSA", "O2S":"OSB", "O3S":"OSC", "C2N":"C7 ",
"O2N":"O7 ", "CME":"C8 ", }
for at in set([ i for v in list(dict_atoms.values()) for i in v]):
    dict_atname[at] = at

residues = [ [] for i in range(nres)]
atindices = [ [] for i in range(natoms)]
struc = 1
for l in pdb:
    if l.startswith("ENDMDL"):
        printmodel(atindices, residues, nres, struc)
        struc+=1
        residues = [ [] for i in range(nres)]
        atindices = [ [] for i in range(natoms)]
    if not l.startswith("HETATM") and not l.startswith("ATOM"):
        continue
    if int(l[8:11]) not in list(dict_at2res.keys()) or l[17:20] == "ROH":
        continue
    atname, res, resname  = dict_atname[l[13:16]], l[17:20], dict_names[ l[17:20] ]
    resindex = dict_at2res[int(l[8:11])]
    l2 = "%s%4s%4s  %4d%s" % (l[:11], atname, resname, resindex+1, l[26:-1])
    residues[ resindex ].append(l2)
    try:
        atindices[ resindex ].append(dict_atoms[resname].index(atname))
    except:
        print((resindex, resname, atname), file=sys.stderr)
        print(atindices, file=sys.stderr)
        print(dict_atoms[resname], file=sys.stderr)
        raise

if not l.startswith("ENDMDL"):
    printmodel(atindices, residues, nres, struc)
