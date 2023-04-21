#!/usr/bin/env python
from phaser import *
import argparse
from Bio.PDB import PDBParser

parser = argparse.ArgumentParser()
parser.add_argument("-m", "--mtz", metavar="mtz", help="MTZ file of interest", type=str)
parser.add_argument("-p", "--pdb1", metavar="pdb1", help="PDB file of interest 1", type=str)
parser.add_argument("-q", "--pdb2", metavar="pdb2", help="PDB file of interest 2", type=str)
args = parser.parse_args()

#p = PDBParser()
#pdb1 = p.get_structure("pdb1", args.pdb1)
#pdb2 = p.get_structure("pdb2", pdb2.pdb2)

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
    i.addCOMP_PROT_MW_NUM(28853,1)
    i.addCOMP_PROT_MW_NUM(17522,1)
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
