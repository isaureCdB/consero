# Example for heparin (IDS-SGN)*3

d=`pwd`
# TODO: give full path to /scripts directory
export FRAG=$d/scripts

LANG=en_US
set -u -e


# If your input is a list of PDB files:
if false;then
##################################################
echo "Create a multi-pdb from AD3 output"
##################################################
cd $d/Docking_solutions
for dp in dp3 dp3a; do
    cd HE_$dp
    j=1
    cat /dev/null >  ../../$dp.pdb
    for i in `ls *.pdb`; do
        echo "MODEL $j" >> ../../$dp.pdb
        grep ATOM $i >> ../../$dp.pdb
        echo "ENDMDL" >> ../../$dp.pdb
        j=$(($j+1))
    done
    cd ../
done
cd ../
fi

# If your input is a numpy array of coordinates:
cd $d/Docking_solutions
for dp in dp3 dp3a; do
    $FRAG/npy2pdb.py HE_$dp.npy template_HE_$dp.pdb > ../$dp.pdb
done
ls ../dp3*-aa.pdb > ../pdbfiles.list
cd $d

##################################################
echo "convert into all-atom / reduced formats"
##################################################
# g1-3.pdb = residues 1 to 3 in the bound ligand.
# g1-3r.pdb : r="reduced"
# TODO : arrange wich one of dp3/dp3a corresponf to g1-3.pdb or g2-4.pdb
#
# convert into same format as bound form (atom names, res names, atom order...)
# args: [docking poses] [atom indices] [template PDB]
# see in ./reduce-gl_indices.py for details on the inputs
$FRAG/reduce-gl_indices.py dp3a.pdb dp3a-aa.ind boundfrag/g2-4-aa.pdb > dp3a-aa.pdb
$FRAG/reduce-gl_indices.py dp3.pdb dp3-aa.ind  boundfrag/g1-3-aa.pdb > dp3-aa.pdb

# convert into "reduced" form : only S, O in cycles and C of COO atoms
$FRAG/reduce-gl_indices.py dp3a.pdb dp3ar.ind boundfrag/g2-4r.pdb > dp3ar.pdb
$FRAG/reduce-gl_indices.py dp3.pdb dp3r.ind  boundfrag/g1-3r.pdb > dp3r.pdb
#$FRAG/reduce-gl_indices.py dp3a-1.pdb dp3a-cg.ind boundfrag/g2-4r.pdb > dp3ar-1.pdb

##################################################
echo "convert PDB into numpy arrays of coordinates"
##################################################
# $dp.preatoms:
#   indices (from 1) of the overlapping atoms with next fragment
#   (= atoms in the 2 last residues)
#
# $dp.postatoms:
#   indices (from 1) of the overlapping atoms with previous fragment
#  (= atoms in the 2 first residues)
#
# convert $dp-cg.ind into $dp.preatoms
$FRAG/ind2prepostatom.sh dp3ar.ind dp3ar
$FRAG/ind2prepostatom.sh dp3r.ind dp3r
$FRAG/ind2prepostatom.sh dp3a-aa.ind dp3a-aa
$FRAG/ind2prepostatom.sh dp3-aa.ind dp3-aa

# preatoms.npy = coordinates of the 2 last residues for all poses
# postatoms.npy = coordinates of the 2 first residues for all poses
# format of pre/postatoms: numpy array of shape [nb poses; nb atoms; 3 coordinates]
for dp in dp3 dp3a; do
    $FRAG/pdb2npy.py $dp-aa.pdb --outp $dp-aa.npy
    $FRAG/pdb2npy.py $dp\r.pdb --outp $dp\r.npy
    $FRAG/select-npy.py $dp\r.npy $dp-preatoms.npy --atoms `cat $dp\r.preatoms`
    $FRAG/select-npy.py $dp\r.npy $dp-postatoms.npy --atoms `cat $dp\r.postatoms`
done

link(){
  dp=$1
  i=$2
  for atom in preatoms postatoms; do
    ln -s $dp-aa.npy frag$i-aa.npy
    ln -s $dp-aa.$atom frag$i-aa.$atom
    ln -s $dp\r.$atom frag$i\r.$atom
    ln -s $dp-$atom.npy frag$i-$atom.npy
  done
}

set -u +e

for i in 1 3 5; do
  link dp3 $i
done

for i in 2 4 6; do
  link dp3a $i
done
set -u -e

######################################################################
echo "compute the RMSD of each pose toward the bound structure"
######################################################################
$FRAG/rmsdnpy.py dp3r.npy --ref boundfrag/g1-3r.pdb > frag1.rmsd
$FRAG/rmsdnpy.py dp3r.npy --ref boundfrag/g3-5r.pdb > frag3.rmsd
#rmsdnpy.py dp3r.npy boundfrag/g5-7r.pdb > frag5.rmsd

$FRAG/rmsdnpy.py dp3ar.npy --ref boundfrag/g2-4r.pdb > frag2.rmsd
$FRAG/rmsdnpy.py dp3ar.npy --ref boundfrag/g4-6r.pdb > frag4.rmsd
