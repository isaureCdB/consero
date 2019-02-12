d=`pwd`
export FRAG=$d/scripts

# pdbfiles.list = list of PDB files of docking poses
pdbfiles=$1
#clustering cutoff in Angström (advised: ~0.5)
cutoff=$2

mkdir monomers_library
cd monomers_library
set -u -e
$FRAG/extract-monomers.py ../$pdbfiles > monomers.list

for motif in `cat monomers.list` ;do
    #remove very similar fragments after fitting
    $FRAG/deredundant-npy.py $motif.npy 0.3 --outp $motif-dr.npy --fit --center
    #cluster at $cutoff Angström
    $FRAG/npy_pairwise_rmsd_cluster.py $motif-dr.npy $cutoff
    rm $motif.npy $motif-dr.npy
    ln -s $motif-dr-clust$cutoff.npy $motif.npy
done

cd ../
