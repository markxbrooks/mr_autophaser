#!/usr/bin/env ccp4-python
"""
mr_autophaser.py
(c) 2023 Mark Brooks
Takes 1 or more input pdb files, performs molecular replacement 
using a given MTZ file
"""
import os,sys
import re
from phaser import *
from StringIO import *
from optparse import OptionParser
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from Bio.SeqUtils.ProtParam import ProteinAnalysis


def main(options, args):
    CHAIN = options.CHAIN
    MTZIN = options.MTZIN
    pdb_list = {options.PDBIN1: options.NUMBER1, 
            options.PDBIN2: options.NUMBER2, 
            options.PDBIN3: options.NUMBER3}
    i = InputMR_DAT()
    i.setHKLI(MTZIN)
    i.setMUTE(True)
    r = runMR_DAT(i)
    # print(r.logfile())
    print("Running MR Auto-Phaser...")
    if r.Success():
        i = InputMR_AUTO()
        i.setREFL_DATA(r.getREFL_DATA())
        rootname = ""
        for pdb, num in pdb_list.items():
            if pdb:
                fileroot = re.sub(".pdb", "", pdb, re.IGNORECASE)
                if rootname == "":
                    rootname = fileroot    
                else:    
                    rootname = fileroot + "_" + rootname
                i.addENSE_PDB_ID(pdb, pdb,1.0)
                pdbparser = PDBParser(QUIET=True)
                structure = pdbparser.get_structure(pdb, pdb)
                chains = {chain.id:seq1(''.join(residue.resname for residue in chain)) for chain in structure.get_chains()}
                # print("Chains: " + str(chains))
                query_chain = chains[CHAIN]
                # print(query_chain)
                prot_param = ProteinAnalysis(query_chain)
                mwt = prot_param.molecular_weight()
                print("Molecular weight of " + pdb + " Chain '" + CHAIN + "': %0.2f Da" % mwt)
                i.addCOMP_PROT_MW_NUM(mwt,int(num))
                i.addSEAR_ENSE_NUM(pdb,int(num))
        i.setROOT(rootname)
        i.setMUTE(True)
        del(r)
        r = runMR_AUTO(i)
        if r.Success():
            if r.foundSolutions() :
                print("Phaser has found " + str(r.numSolutions()) + " MR solutions")
                print("Top LLG = %f" % r.getTopLLG())
                print("Top TFZ = %f" % r.getTopTFZ())
                print("Top PDB file = %s" % r.getTopPdbFile())
            else:
                print("Phaser has not found any MR solutions")
        else:
            print("Job exit status FAILURE")
            print(r.ErrorName(), "ERROR: ", r.ErrorMessage())
    else:
        print("Job exit status FAILURE")
        print(r.ErrorName(), "ERROR: ", r.ErrorMessage())

if __name__ == '__main__':
    usage = 'This script takes 1 or more input pdb files, performing molecular replacement using a given MTZ file \n./mr_autophaser.py -m beta_blip_P3221.mtz -1 beta.pdb -2 blip.pdb -c " "'
    parser = OptionParser(usage=usage)
    parser.add_option("-m", "--mtz_input", dest="MTZIN",
                      help="MTZ input file", metavar="MTZIN")
    parser.add_option("-1", "--pdb_input1", dest="PDBIN1",
                      help="pdb input file", metavar="PDBIN1")
    parser.add_option("-2", "--pdb_input2", dest="PDBIN2",
                      help="pdb input file", metavar="PDBIN2")
    parser.add_option("-3", "--pdb_input3", dest="PDBIN3",
                      help="pdb input file3", metavar="PDBIN3")
    parser.add_option("-c", "--chain", dest="CHAIN",
                      help="chain name, e.g. 'A','C',' ' Default is 'A'", metavar="CHAIN", default="A")
    parser.add_option("-n", "--num_pdb1", dest="NUMBER1",
                      help="-num_pdb1, default is 1", metavar="NUMBER1", default=1)
    parser.add_option("-o", "--num_pdb2", dest="NUMBER2",
                      help="-num_pdb2, default is 1", metavar="NUMBER2", default=1)
    parser.add_option("-p", "--num_pdb3", dest="NUMBER3",
                      help="-num_pdb3, default is 1", metavar="NUMBER3", default=1)
    (options, args) = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        print "No argument given!"
        parser.print_help()
    else:    
        main(options, args)
