#!/usr/bin/env ccp4-python
"""
mr_autophaser.py
(c) 2023 Mark Brooks
Takes 1 or more input pdb files, performs molecular replacement 
using a given MTZ file
ccp4-python ./mr_autophaser.py -m test/beta_blip.mtz -1 test/beta.pdb -2 test/blip.pdb
"""
import sys
from pathlib import Path
from time import strftime
import logging
import argparse
from phaser import *
from Bio.PDB import PDBParser, PPBuilder
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from utils import print_header

__version__ = "0.1.0"
__program__ = "mr_autophaser"
__author__ = "Mark Brooks, based on phaser tutorial scripts"
__synopsis__ = """This script takes 1 or more input pdb files, performing molecular
replacement using a given MTZ file 
example usage:
ccp4-python ./mr_autophaser.py -m test/beta_blip.mtz -1 test/beta.pdb -2 test/blip.pdb
"""

def parse_args():
    """
    @synopsis parse arguments
    @return args
    """
    parser = argparse.ArgumentParser(usage=__synopsis__)
    if len(sys.argv[1:]) == 0:
        print("No argument given!")
        parser.print_help()
    parser.add_argument(
        "-m", "--mtz_input", dest="mtzin", help="MTZ input file", metavar="mtzin"
    )
    parser.add_argument(
        "-1", "--pdb_input1", dest="PDBIN1", help="pdb input file", metavar="PDBIN1"
    )
    parser.add_argument(
        "-2", "--pdb_input2", dest="PDBIN2", help="pdb input file", metavar="PDBIN2"
    )
    parser.add_argument(
        "-3", "--pdb_input3", dest="PDBIN3", help="pdb input file3", metavar="PDBIN3"
    )
    parser.add_argument(
        "-c",
        "--chain",
        dest="chain",
        help="chain name, e.g. 'A','C',' ' Default is 'None'",
        metavar="chain",
        default="",
    )
    parser.add_argument(
        "-n",
        "--num_pdb1",
        dest="NUMBER1",
        help="-num_pdb1, default is 1",
        metavar="NUMBER1",
        default=1,
    )
    parser.add_argument(
        "-o",
        "--num_pdb2",
        dest="NUMBER2",
        help="-num_pdb2, default is 1",
        metavar="NUMBER2",
        default=1,
    )
    parser.add_argument(
        "-p",
        "--num_pdb3",
        dest="NUMBER3",
        help="-num_pdb3, default is 1",
        metavar="NUMBER3",
        default=1,
    )
    args = parser.parse_args()
    return args


def main(args):
    """
    Set up logging
    Parse arguments
    Perform molecular replacement
    @param args
    """
    print_header(__program__, __author__, __version__, __synopsis__)
    timestamp = strftime("%H_%M_%m_%d_%Y")
    log_filename = f"logs/{__program__}_{timestamp}.log"
    print(f"Saving logs to: {log_filename}")
    logging.basicConfig(
        filename=log_filename,
        filemode="w",
        format='%(asctime)s,%(msecs)d %(levelname)-8s [%(pathname)s:%(lineno)d in function %(funcName)s] %(message)s',
        level=logging.DEBUG,
    )
    chain = args.chain
    mtzin = args.mtzin
    pdb_list = {
        args.PDBIN1: args.NUMBER1,
        args.PDBIN2: args.NUMBER2,
        args.PDBIN3: args.NUMBER3,
    }
    mr_auto(mtzin, pdb_list, chain)

def mr_auto(mtzin, pdb_list, chain=None):
    """
    @synopsis perform molecular replacement 
    @param pdb_list dictionary of pdb files and the number of copies of each expected in the unit cell 
    @param chain string chain of interest
    @return True on success
    """
    phaser_input = InputMR_DAT()
    phaser_input.setHKLI(mtzin)
    phaser_input.setMUTE(True)
    phaser_result = runMR_DAT(phaser_input)
    logging.info(phaser_result.logfile())
    print("Running MR Auto-Phaser...")
    if phaser_result.Success():
        phaser_input = InputMR_AUTO()
        phaser_input.setREFL_DATA(phaser_result.getREFL_DATA())
        rootname = ""
        for pdb, num in pdb_list.items():
            if pdb:
                fileroot = Path(pdb).stem
                if rootname == "":
                    rootname = fileroot
                else:
                    rootname = fileroot + "_" + rootname
                phaser_input.addENSE_PDB_ID(pdb, pdb, 1.0)
                pdbparser = PDBParser(QUIET=True)
                if chain:
                    structure = pdbparser.get_structure(pdb, pdb)
                    chains = {
                        chain.id: seq1("".join(residue.resname for residue in chain))
                        for chain in structure.get_chains()
                    }
                    query_chain = chains[chain]
                    # print(query_chain)
                    prot_param = ProteinAnalysis(query_chain)
                    mwt = prot_param.molecular_weight()
                    print(f">{fileroot}")
                    print((query_chain))
                    print(f"Molecular weight of {pdb} Chain {chain}: {mwt:0.2f} Da")
                else:
                    structure = pdbparser.get_structure(pdb, pdb)
                    pdbb = PPBuilder()
                    for peptide in pdbb.build_peptides(structure):
                        query_chain = peptide.get_sequence().__str__()
                    analysis = ProteinAnalysis(str(query_chain))
                    mwt = analysis.molecular_weight()
                    print(f">{fileroot}")
                    print(query_chain)
                    print(f"Molecular weight of {pdb} {mwt:0.2f} Da")
                phaser_input.addCOMP_PROT_MW_NUM(mwt, int(num))
                phaser_input.addSEAR_ENSE_NUM(pdb, int(num))
        phaser_input.setROOT(rootname)
        phaser_input.setMUTE(True)
        del phaser_result
        phaser_result = runMR_AUTO(phaser_input)
        if phaser_result.Success():
            if phaser_result.foundSolutions():
                print(
                    f"Phaser has found {str(phaser_result.numSolutions())} MR solutions"
                )
                print(f"Top LLG = {phaser_result.getTopLLG():0.2f}")
                print(f"Top TFZ = {phaser_result.getTopTFZ():0.2f}")
                print(f"Top PDB file = {phaser_result.getTopPdbFile()}")
                return True
            else:
                print("Phaser has not found any MR solutions")
                return False
        else:
            print("Job exit status FAILURE")
            print(phaser_result.ErrorName(), "ERROR: ", phaser_result.ErrorMessage())
            return False
    else:
        print("Job exit status FAILURE")
        print(phaser_result.ErrorName(), "ERROR: ", phaser_result.ErrorMessage())
        return False    
    logging.info(phaser_result.logfile())

if __name__ == "__main__":
    args = parse_args()
    main(args)
