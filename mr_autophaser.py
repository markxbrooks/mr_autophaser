#!/usr/bin/env python
from pathlib import Path
from phaser import *
import argparse
from Bio.PDB import PDBParser, PPBuilder
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

parser = argparse.ArgumentParser()
parser.add_argument("-m", "--mtz", metavar="mtz", help="MTZ file of interest", type=str)
parser.add_argument("-p", "--pdb1", metavar="pdb1", help="PDB file of interest 1", type=str)
parser.add_argument("-q", "--pdb2", metavar="pdb2", help="PDB file of interest 2", type=str)
args = parser.parse_args()

p1_stem = Path(args.pdb1).stem
p2_stem = Path(args.pdb2).stem

p = PDBParser(QUIET=True)
pdb1_structure = p.get_structure("pdb1", args.pdb1)
ppb1 = PPBuilder()
for pp1 in ppb1.build_peptides(pdb1_structure):
    record1_str = pp1.get_sequence().__str__()
    
pdb2_structure = p.get_structure("pdb2", args.pdb2)
ppb2 = PPBuilder()
for pp2 in ppb2.build_peptides(pdb2_structure):
    record2_str = pp2.get_sequence().__str__()
    
print(f">{p1_stem}")
print(record1_str)
print(f">{p2_stem}")
print(record2_str)

p1_analysis = ProteinAnalysis(str(record1_str))
p1_mr = p1_analysis.molecular_weight()

p2_analysis = ProteinAnalysis(str(record2_str))
p2_mr = p2_analysis.molecular_weight()

print(f"{p1_stem} pdb1 mr: {p1_mr:.1f} Da")
print(f"{p2_stem} pdb2 mr: {p2_mr:.1f} Da")

i = InputMR_DAT()
i.setHKLI(args.mtz)
i.setLABI_F_SIGF("Fobs","Sigma")
i.setHIRES(6.0)
i.setMUTE(True)
r = runMR_DAT(i)
if r.Success():
    i = InputMR_AUTO()
    i.setSPAC_HALL(r.getSpaceGroupHall())
    i.setCELL6(r.getUnitCell())
    i.setREFL_F_SIGF(r.getMiller(),r.getFobs(),r.getSigFobs())
    i.setROOT("beta_blip_auto")
    i.addENSE_PDB_ID("pdb1",args.pdb1,1.0)
    i.addENSE_PDB_ID("pdb2",args.pdb2,1.0)
    i.addCOMP_PROT_MW_NUM(p1_mr,1)
    i.addCOMP_PROT_MW_NUM(p2_mr,1)
    i.addSEAR_ENSE_NUM("pdb1",1)
    i.addSEAR_ENSE_NUM("pdb2",1)
    i.setMUTE(True)
    del(r)
    r = runMR_AUTO(i)
    if r.Success():
        if r.foundSolutions() :
            print("Phaser has found MR solutions")
            print("Top LLG = %f" % r.getTopLLG())
            print("Top PDB file = %s" % r.getTopPdbFile())
        else:
            print("Phaser has not found any MR solutions")
    else:
        print("Job exit status FAILURE")
        print(r.ErrorName(), "ERROR :", r.ErrorMessage())
else:
    print("Job exit status FAILURE")
    print(r.ErrorName(), "ERROR :", r.ErrorMessage())
