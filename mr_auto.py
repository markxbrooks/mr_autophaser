#!/usr/bin/env ccp4-python
"""
mr_autophaser.py
(c) 2023 Mark Brooks
Takes 1 or more input pdb files, performs molecular replacement 
    using a given MTZ file
    ccp4-python ./mr_autophaser.py -m test/beta_blip.mtz -1 test/beta.pdb -2 test/blip.pdb
"""
from pathlib import Path
import logging
from phaser import *
from Bio.PDB import PDBParser, PPBuilder
from Bio.SeqUtils.ProtParam import ProteinAnalysis

__version__ = "0.1.0"
__program__ = "mr_autophaser"
__author__ = "Mark Brooks, based on phaser tutorial scripts"
__synopsis__ = """This script takes 1 or more input pdb files, performing molecular
    replacement using a given MTZ file 
    example usage:
    ccp4-python ./mr_autophaser.py -m test/beta_blip.mtz -1 test/beta.pdb -2 test/blip.pdb
"""

def mr_auto(mtzin, pdb_list, chain=None, rootname=""):
    """
    @synopsis perform molecular replacement 
    @param pdb_list dictionary of pdb files and the number of 
        copies of each expected in the unit cell 
    @param chain string chain of interest
    @return True on success
    """
    phaser_input = InputMR_DAT()
    phaser_input.setHKLI(mtzin)
    phaser_input.setMUTE(True)
    phaser_run = runMR_DAT(phaser_input)
    logging.info(phaser_run.logfile())
    print("Running MR Auto-Phaser...")
    if phaser_run.Success():
        phaser_input = InputMR_AUTO()
        phaser_input.setREFL_DATA(phaser_run.getREFL_DATA())
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
                    sequence = ""
                    chains = {
                        chain.id: sequence("".join(residue.resname for residue in chain))
                        for chain in structure.get_chains()
                    }
                    query_chain = chains[chain]
                    prot_param = ProteinAnalysis(query_chain)
                    mwt = prot_param.molecular_weight()
                    print(f">{fileroot}")
                    print((query_chain))
                    print(f"Molecular weight of {pdb} Chain {chain}: {mwt:0.2f} Da")
                else:
                    structure = pdbparser.get_structure(pdb, pdb)
                    pdbb = PPBuilder()
                    for peptide in pdbb.build_peptides(structure):
                        query_chain = peptide.get_sequence()
                    analysis = ProteinAnalysis(str(query_chain))
                    mwt = analysis.molecular_weight()
                    print(f">{fileroot}")
                    print(query_chain)
                    print(f"Molecular weight of {pdb} {mwt:0.2f} Da")
                phaser_input.addCOMP_PROT_MW_NUM(mwt, int(num))
                phaser_input.addSEAR_ENSE_NUM(pdb, int(num))
        phaser_input.setROOT(rootname)
        phaser_input.setMUTE(True)
        del phaser_run
        phaser_run = runMR_AUTO(phaser_input)
        if phaser_run.Success():
            if phaser_run.foundSolutions():
                print(
                    f"Phaser has found {str(phaser_run.numSolutions())} MR solutions"
                )
                print(f"Top LLG = {phaser_run.getTopLLG():0.2f}")
                print(f"Top TFZ = {phaser_run.getTopTFZ():0.2f}")
                print(f"Top PDB file = {phaser_run.getTopPdbFile()}")
                logging.info(phaser_run.logfile())
                return True
            else:
                print("Phaser has not found any MR solutions")
                logging.info(phaser_run.logfile())
                return False
        else:
            print("Job exit status FAILURE")
            print(phaser_run.ErrorName(), "ERROR: ", phaser_run.ErrorMessage())
            logging.info(phaser_run.logfile())
            return False
    else:
        print("Job exit status FAILURE")
        print(phaser_run.ErrorName(), "ERROR: ", phaser_run.ErrorMessage())
        logging.info(phaser_run.logfile())
        return False
