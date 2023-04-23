#!/usr/bin/env ccp4-python
"""
mr_autophaser.py
(c) 2023 Mark Brooks
Takes 1 or more input pdb files, performs molecular replacement 
    using a given MTZ file
    ccp4-python ./mr_autophaser.py -m test/beta_blip.mtz -1 test/beta.pdb -2 test/blip.pdb
"""
import sys
from time import strftime
import logging
import argparse
from utils import print_header
from mr_auto import mr_auto

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
        sys.exit(-1)
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
        format='%(asctime)s,%(msecs)d %(levelname)-8s \
            [%(pathname)s:%(lineno)d in function %(funcName)s] %(message)s',
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

if __name__ == "__main__":
    arguments = parse_args()
    main(arguments)
