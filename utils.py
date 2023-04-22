#!/bin/env python3.8
"""
Synopsis: Utility functions 
"""

__program__ =  "mxpipe_utils"
__version__ = 0.1
__author__ = "Mark Brooks"
__synopsis__ = "Utility functions for mxpipe"


def print_header(program, author, version, synopsis):
    """Prints header for the file of interest"""
    header=f"""
###########################################################################################

    {program} {version} {author} 
    {synopsis}
###########################################################################################
"""
    print(header)
