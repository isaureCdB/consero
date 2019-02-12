#!/usr/bin/env python2
# Copyright (C) Isaure Chauvot de Beauchene & Sjoerd J. de Vries (TUM)

import sys, numpy as np, weakref
import json, itertools, bisect, argparse

def npy3to2(npy):
    if len(npy.shape) == 3:
        assert npy.shape[2] == 3, npy.shape
        npy = npy.reshape(npy.shape[0], 3*npy.shape[1])
    else:
        assert len(npy.shape) == 2 and npy.shape[1]%3 == 0
    return npy

MINCHILD = 100
MAXCHUNK = 4000000

###########
#$FRAG/assemble.py $nfrag $rmsd $maxstruc $a $b $c > $graph.json
parser =argparse.ArgumentParser(description=__doc__,
                        formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--nfrag', help="number of fragments to assemble")
parser.add_argument('--rmsd', help="maximal overlap RMSD for consecutive fragments")
#parser.add_argument('--top', help="number of top-ranked chain to consider", action="store_true")
parser.add_argument('--maxstruct', help="maximal number of poses per fragment to consider")
parser.add_argument('--preatoms', nargs='+', help="frag1.preatoms frag2.preatoms...")
parser.add_argument('--postatoms', nargs='+', help="frag1.postatoms frag2.postatoms...")
parser.add_argument('--selections', nargs='+', help="frag1.sel frag2.sel...", default=None)

args = parser.parse_args()
############

nfrags = int(args.nfrag)
max_rmsd = float(args.rmsd)**2
maxstruc = int(args.maxstruct)

preatoms = [npy3to2(np.load(f)) for f in args.preatoms]
postatoms = [npy3to2(np.load(f)) for f in args.postatoms]

assert len(postatoms) == len(preatoms)

print >> sys.stderr, "NFRAGS", nfrags
assert nfrags >= 2

CLUSTERING = [10, 8, 6, 5, 4, 3.5, 3, 2.5, 2, 1.7, 1.5, 1.0, 0.5, 0.1, 0]
#last CLUSTERING must be 0, second-to-last is deredundant criterion
CLUST_MARGIN = 2 #2 is the theoretical worse case; change to 9999 to disable all effects of clustering

selections = [[] for n in range(nfrags)]
if args.selections is not None:
  selections = [np.array(sorted([int(l) for l in open(f) if len(l.strip())])) for f in args.selections]
  assert len(selections) == len(postatoms)

for arr in (postatoms, preatoms):
  for anr, a in enumerate(arr):
    ncoor = a.shape[1] / 3
    arr[anr] = a.reshape(a.shape[0], ncoor, 3)

#take top
if maxstruc > 0:
  for arr in (postatoms, preatoms):
    for anr, a in enumerate(arr):
      arr[anr] = arr[anr][:maxstruc]
  for selnr, sel in enumerate(selections):
    if not len(sel): continue
    pos = bisect.bisect_right(sel, maxstruc)
    selections[selnr] = sel[:pos]
    assert len(selections[selnr])

nstruc = []
for n in range(nfrags):
  a1 = postatoms[n]
  a2 = preatoms[n]
  assert a1.shape[0] == a2.shape[0], (n, a1.shape, a2.shape)
  nstruc.append(a1.shape[0])

for n in range(1, nfrags):
  a1 = postatoms[n-1]
  a2 = preatoms[n]
  assert a1.shape[1] == a2.shape[1], (n, a1.shape, n+1, a2.shape)

for selnr, sel in enumerate(selections):
  for conf in sel:
    assert conf > 0 and conf <= nstruc[selnr], (conf, nstruc[selnr])

ranks = [np.arange(s)+1 for s in nstruc]

for n in range(nfrags):
  if len(selections[n]):
    postatoms[n] = postatoms[n][selections[n]-1]
    preatoms[n] = preatoms[n][selections[n]-1]
    ranks[n] = selections[n]
    nstruc[n] = len(selections[n])

MAX_SIZEKEY = 100000000000
MAX_CLUSTERING = len(CLUSTERING) - 1
class Cluster(object):
  __slots__ = ("clustid", "_splittable", "clusterlevel", "coors", "ranks", "all_ranks", "children", "nodes", "totstruc", "parent", "connections", \
    "back_connections", "_splitting", "_checking_delete", "conlevel")
  def __init__(self, clustid, clusterlevel, coors, ranks):
    self.clustid = clustid
    self._splittable = True
    if coors is not None and len(coors) == 1: #singleton
      self.clusterlevel = MAX_CLUSTERING
    else:
      assert clusterlevel < MAX_CLUSTERING
      self.clusterlevel = clusterlevel #contains the clusterlevel of the cluster itself, not of the cluster children!
    if self.clusterlevel == MAX_CLUSTERING:
      self._splittable = False
    self.coors = coors
    self.ranks = ranks
    self.all_ranks = set(ranks)
    self.children = []
    self.nodes = 1
    if coors is not None:
      self.totstruc = coors.shape[0]
    self.parent = None
    self.connections = []
    self.back_connections = []
    self._splitting = False
    self._checking_delete = False
    self.conlevel = self.clusterlevel
    if self.clusterlevel is None:
      self.conlevel = -1
    if coors is not None and clusterlevel == MAX_CLUSTERING - 1: #deredundant level
      r = self
      for cnr in range(len(self.coors)):
        c = Cluster(self.clustid + (cnr,), MAX_CLUSTERING, coors[cnr:cnr+1], ranks[cnr:cnr+1])
        c.parent = r
        self.children.append(c)
      self.nodes = len(self.children)
      self._splittable = False

  def _cluster(self, clusterlevel):
    assert not self.children

    c = self.coors
    #clus: coordinates of the first structure of each cluster
    clus = c[:1]
    #clus_indices: the coors indices of the structures of each cluster
    clus_indices = [[0]]
    chunksize = 20

    radius = CLUSTERING[clusterlevel]
    max_sd = radius * radius * c.shape[1]

    #This variable keeps track, for every structure in the chunk, into which new cluster it is sorted
    which_new_clust = np.zeros(chunksize, dtype=int)

    clustid = self.clustid
    if clustid is None: clustid = ()
    for n in range(1, len(c), chunksize):
      chunk = c[n:n+chunksize]

      #intra-chunk msd
      d = chunk[:, np.newaxis, :, :] - chunk[np.newaxis, :, :, :]
      intra_msd = np.einsum("...ijk,...ijk->...i", d,d)

      #chunk-cluster msd
      d = chunk[:, np.newaxis, :, :] - clus[np.newaxis, :, :, :]
      inter_msd = np.einsum("...ijk,...ijk->...i", d,d)

      for nn in range(len(chunk)):
        sort_clust = None
        close_intra_clusts = (intra_msd[nn] < max_sd).nonzero()[0]
        intra_new_clusts = [which_new_clust[k] for k in close_intra_clusts if k < nn and which_new_clust[k] != -1]
        if len(intra_new_clusts):
          sort_clust = min(intra_new_clusts)
        close_inter_clusts = (inter_msd[nn] < max_sd).nonzero()[0]
        if len(close_inter_clusts):
          sort_clust2 = min(close_inter_clusts)
          if sort_clust is None or sort_clust > sort_clust2:
            sort_clust = sort_clust2
        if sort_clust is None:
          #new cluster
          which_new_clust[nn] = len(clus)
          clus = np.append(clus, chunk[nn][np.newaxis,:,:], axis=0)
          clus_indices.append([n+nn])
        else:
          clus_indices[sort_clust].append(n+nn)
          which_new_clust[nn] = -1

    indices = [i[0] for i in clus_indices]

    #Re-cluster to the lowest RMSD cluster
    clus_indices = [i[:1] for i in clus_indices]
    for n in range(0, len(c), chunksize):
      chunk = c[n:n+chunksize]
      d = chunk[:, np.newaxis, :, :] - clus[np.newaxis, :, :, :]
      inter_msd = np.einsum("...ijk,...ijk->...i", d,d)
      sort_clusts = np.argmin(inter_msd, axis=1)
      for nn in range(len(chunk)):
        if (n+nn) in indices: continue
        sort_clust = sort_clusts[nn]
        clus_indices[sort_clust].append(n+nn)


    for cnr,c in enumerate(clus_indices):
      ind = clus_indices[cnr]
      coors = self.coors[ind]
      ranks = self.ranks[ind]
      c = Cluster(clustid+(cnr+1,), clusterlevel, coors, ranks)
      self.children.append(c)
    self.nodes = len(self.children)
    self.totstruc = sum([c.totstruc for c in self.children])
    return clus, indices
  def cluster(self, clusterlevel):
    assert clusterlevel < MAX_CLUSTERING
    clus, indices = self._cluster(clusterlevel)
    self.coors = clus
    self.ranks = self.ranks[indices]
    r = self
    for c in self.children:
      c.parent = r
    for c in self.children:
      if c._splittable:
        break
    else:
      self._splittable = False
  def dissolve(self, clusterlevel):
    self.ranks = np.concatenate([c.ranks for c in self.children])
    newchildren = []
    coors = []
    while len(self.children):
      child = self.children.pop()
      child.cluster(clusterlevel)
      coors.append(child.coors)
      newchildren += child.children
      self.totstruc = child.totstruc
    self.coors = np.concatenate(coors, axis=0)
    self.children = newchildren
    r = self
    for c in self.children:
      c.parent = r
    self.nodes = len(newchildren)
    for c in self.children:
      if c._splittable:
        break
    else:
      self._splittable = False
  def reorganize(self):
    for c in self.children:
      c.reorganize()
    if not len(self.children): return
    if len(self.children) >= MINCHILD: return
    oldchildren = [c for c in self.children if len(c.children)]
    if not len(oldchildren): return
    while len(self.children) < MINCHILD and len(oldchildren):
      child = oldchildren.pop(0)
      self.children.remove(child)
      self.children += child.children
      oldchildren += [c for c in child.children if len(c.children)]
    r = self
    coors = [c.coors[0] for c in self.children]
    self.coors = np.array(coors)
    for c in self.children:
      c.parent = r
    self.nodes = sum([c.nodes for c in self.children])
    for c in self.children:
      if c._splittable:
        break
    else:
      self._splittable = False

  def split(self):
    if not self.children:
      assert self.clusterlevel is not None
      if self.clusterlevel == MAX_CLUSTERING: return False
      self.clusterlevel += 1
      if self.clusterlevel == MAX_CLUSTERING:
        self._splittable = False
        return False
      self.cluster(self.clusterlevel)
      r = self
      while len(self.children) == 1:
        child = self.children[0]
        self.clusterlevel = child.clusterlevel
        if self.clusterlevel >= MAX_CLUSTERING - 1:
          self._splittable = False
          break
        self.clusterlevel += 1
        child.cluster(self.clusterlevel)
        children = child.children
        for c in children: c.parent = r
        self.coors = child.coors
        self.ranks = child.ranks
        self.children = children
        self.nodes = len(children)
        self.totstruc = sum([c.totstruc for c in children])

      self.parent.add_nodes(self.nodes - 1)
      for c in self.children:
        if c._splittable:
          break
      else:
        self._splittable = False
      return True
    else:
      self._splitting = True
      oldnodes = self.nodes
      ok = False
      for c in self.children:
        c_splittable = c._splittable
        has_split = c.split()
        if has_split: ok = True
      newnodes = self.nodes
      self._splitting = False
      if self.parent is not None and newnodes > oldnodes:
        self.parent.add_nodes(newnodes-oldnodes)
      for c in self.children:
        if c._splittable:
          break
      else:
        self._splittable = False
      return ok

  def add_nodes(self, nodes):
    self.nodes += nodes
    if self.parent is not None and not self._splitting:
      self.parent.add_nodes(nodes)

  def check_deletion(self):
    if self._checking_delete: return
    self._checking_delete = True
    while 1:
      has_c1, has_c2 = len(self.back_connections), len(self.connections)
      if has_c1 and has_c2: break
      frag0 = self.clustid[0]
      if frag0 > 1000:
        pos = 2 * (frag0-1001) + 1
      else:
        pos = 2 * (frag0-1)
      if pos == 0 and has_c2: break
      if pos == len(clusters) - 1 and has_c1: break
      if self not in clusters[pos]: break
      clusters[pos].remove(self)
      for o in list(self.back_connections):
        o.connections.remove(self)
        o.check_deletion()
      for o in list(self.connections):
        o.back_connections.remove(self)
        o.check_deletion()
      break
    self._checking_delete = False


  def decompose(self, fwd):
    c1 = self.coors
    if fwd:
      if not self.connections: return
      others = list(self.connections)
      c2 = np.array([c.coors[0] for c in self.connections])
    else:
      if not self.back_connections: return
      others = list(self.back_connections)
      c2 = np.array([c.coors[0] for c in self.back_connections])

    chunksize = MAXCHUNK/len(c1)
    for chunkpos in range(0, len(c2), chunksize):
      c2_chunk = c2[chunkpos:chunkpos+chunksize]
      others_chunk = others[chunkpos:chunkpos+chunksize]
      d = c1[:, np.newaxis, :, :] - c2_chunk[np.newaxis, :, :, :]

      o_max_rmsd = []
      for o in others_chunk:
        if o.clusterlevel is None:
          mx = 10**6
        else:
          mx = CLUSTERING[o.clusterlevel] * CLUST_MARGIN
        o_max_rmsd.append(mx)
      c_max_rmsd = []
      for child in self.children:
        mx = CLUSTERING[child.clusterlevel] * CLUST_MARGIN
        c_max_rmsd.append(mx)

      max_rmsd0 = np.array(c_max_rmsd)[:, np.newaxis] + np.array(o_max_rmsd)[np.newaxis, :]
      max_rmsd2 = max_rmsd0 + max_rmsd
      max_sd = (max_rmsd2**2) * c1.shape[1]

      msd = np.einsum("...ijk,...ijk->...i", d,d)
      if fwd:
        ocon = [o.back_connections for o in others_chunk]
        childcon = [c.connections for c in self.children]
      else:
        ocon = [o.connections for o in others_chunk]
        childcon = [c.back_connections for c in self.children]
      for o in ocon:
        o.remove(self)

      for childnr, onr in zip(*np.where(msd < max_sd)):
        #print childnr, onr, len(self.children), len(ocon)
        ocon[onr].append(self.children[childnr])
        childcon[childnr].append(others_chunk[onr])

      for o in others:
        o.check_deletion()

  def decompose_intra(self, fwd):
    if fwd:
      if not self.connections: return
      others = list(self.connections)
    else:
      if not self.back_connections: return
      others = list(self.back_connections)

    if fwd:
      ocon = [o.back_connections for o in others]
      childcon = [c.connections for c in self.children]
    else:
      ocon = [o.connections for o in others]
      childcon = [c.back_connections for c in self.children]
    for o in ocon:
      o.remove(self)

    for childnr, child in enumerate(self.children):
      cranks = child.all_ranks
      for onr, o in enumerate(others):
        if cranks.intersection(o.all_ranks):
          ocon[onr].append(child)
          childcon[childnr].append(o)

    for o in others:
      o.check_deletion()

  def verify(self):
    if len(self.children): return
    if self.totstruc > 1: return
    cons = [con for con in self.connections if not len(con.children) and con.totstruc == 1]
    if not len(cons): return
    concoors = np.concatenate([con.coors[:1] for con in cons])
    c2 = self.coors[0]
    c1 = concoors
    d = c1 - c2
    msd = np.einsum("...ijk,ijk->...i", d,d)
    msd_low = (msd<(max_rmsd**2*c2.shape[0]))
    for n in range(len(cons)):
      assert msd_low[n]
    return

  def all_children(self):
    if not len(self.children):
      yield self
    else:
      for cc in self.children:
        for v in cc.all_children():
          yield v


#Build cluster tree
clusters = []
for n in range(nfrags):
  for atoms in (preatoms, postatoms):
    i = n + 1
    if atoms is postatoms: i += 1000
    c = Cluster((i,), None, atoms[n], ranks[n])

    if not (n % 2):
      a = atoms[n]
      r = ranks[n]
      for nn in range(len(a)):
        cc = Cluster((i,nn), MAX_CLUSTERING, a[nn:nn+1], np.array(r[nn:nn+1]))
        c.children.append(cc)
      c.nodes = len(a)
      clusters.append([c])
      continue

    clusterlevel = 0
    c.cluster(clusterlevel)
    for clusterlevel in range(1, len(CLUSTERING)-1):
      if len(c.children) >= MINCHILD: break
      c.dissolve(clusterlevel)
    count = 0
    #assert c.clusterlevel is None

    def split_all(c):
      global count
      if not c._splittable: return
      if not len(c.children):
        ok = c.split()
        count += 1
        if not (count % 5000): print >> sys.stderr, n+1, count,  "/", len(atoms[n])
        if not ok: return
      for cc in c.children:
        split_all(cc)

    split_all(c)
    c.reorganize()
    print >> sys.stderr, n+1, nstruc[n], CLUSTERING[clusterlevel], len(c.children), c.nodes
    #assert c.clusterlevel is None
    clusters.append([c])

#Initialize tree connections, intra-fragment, split-down
for n in range(0, 2 * nfrags, 4):
  c1, c2 = clusters[n][0], clusters[n+1][0]
  #print >> sys.stderr,  n, n+1, len(c1.children), len(c2.children)
  for nn in range(len(c1.children)):
    cc1, cc2 = c1.children[nn], c2.children[nn]
    cc1.connections = [cc2]
    cc2.back_connections = [cc1]

#Initialize tree connections, intra-fragment, non-split-down
for n in range(2, 2 * nfrags, 4):
  c1, c2 = clusters[n][0], clusters[n+1][0]
  c1.connections.append(c2)
  c2.back_connections.append(c1)

#Initialize tree connections, inter-fragment, split-down to non-split-down
for n in range(2, 2 * nfrags, 4):
  c1, c2 = clusters[n-1][0], clusters[n][0]
  for cc in c1.children:
    cc.connections = [c2]
    c2.back_connections.append(cc)
  clusters[n-1] = c1.children

#Initialize tree connections, inter-fragment, non-split-down to split-down
for n in range(4, 2 * nfrags, 4):
  c1, c2 = clusters[n-1][0], clusters[n][0]
  for cc in c2.children:
    cc.back_connections = [c1]
    c1.connections.append(cc)
  clusters[n] = c2.children

clusters[0] = clusters[0][0].children
if len(clusters[-2]) > 1:
  clusters[-1] = clusters[-1][0].children

def decompose(clusnr):
  best_conlevel = None
  best = None
  clus = clusters[clusnr]
  for nn in range(len(clus)):
    c = clus[nn]
    if not len(c.children): continue
    if best_conlevel is None or c.conlevel > best_conlevel:
      best_conlevel = c.conlevel
      best = nn
  if best is None: return False
  c = clus[best]
  clusters[clusnr].remove(c)
  if (clusnr % 2):
    c.decompose(fwd=True)
    c.decompose_intra(fwd=False)
  else:
    c.decompose_intra(fwd=True)
    c.decompose(fwd=False)
  for child in c.children:
    conlevel = 0
    for con in itertools.chain(child.back_connections,child.connections):
      v = con.clusterlevel
      if v is not None and v > conlevel: conlevel = v
    child.conlevel = conlevel * child.clusterlevel
  clusters[clusnr] += c.children

  for cc in c.children:
    cc.check_deletion()
  return True

#Decompose tree
step = 0
to_decompose = []
to_decompose0 = list(range(2, len(clusters), 4))
while len(to_decompose0):
  to_decompose.append(to_decompose0.pop(0))
  if not len(to_decompose0): break
  to_decompose.append(to_decompose0.pop(-1))

for clusnr in to_decompose:
  done1, done2 = False, False
  while 1:
    if not (step % 5):
      print >> sys.stderr, [len(c) for c in clusters]
      print >> sys.stderr, [sum([cc.nodes for cc in c]) for c in clusters]
    step += 1
    if not done1:
      ok1 = decompose(clusnr)
      if not ok1:
        if done2: break
        done1 = True
    if not done2:
      ok2 = decompose(clusnr+1)
      if not ok2:
        if done1: break
        done2 = True

print >> sys.stderr, [len(c) for c in clusters[::2]]
print >> sys.stderr, [sum([cc.nodes for cc in c]) for c in clusters[::2]]

#Verification
for c in clusters[1::2]:
  for cc in c:
    cc.verify()

#Sort clusters
for c in clusters:
  c.sort(key=lambda clus: clus.ranks[0])

#write out tree
tree = {"nfrags":nfrags, "max_rmsd": max_rmsd, "clusters": [], "interactions": []}
for cnr in range(1,len(clusters), 2):
  c = clusters[cnr]
  clus = []
  inter = []
  if cnr < len(clusters)-1:
    nex = clusters[cnr+1]
    nexmap = {}
    for vnr, v in enumerate(nex):
      nexmap[v] = vnr
  for ccnr, cc in enumerate(c):
    cclus = {"radius":CLUSTERING[cc.clusterlevel], "ranks": cc.ranks.tolist()}
    clus.append(cclus)
    if cnr < len(clusters)-1:
      for other in cc.connections:
        index = nexmap[other]
        inter.append((ccnr, index))
  inter.sort(key=lambda i: 100000*i[0]+i[1])
  tree["clusters"].append(clus)
  if cnr < len(clusters)-1:
    tree["interactions"].append(inter)

json.dump(tree, sys.stdout)
