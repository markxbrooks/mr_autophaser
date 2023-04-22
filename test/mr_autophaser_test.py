from unittest import TestCase
from mr_autophaser import *

def test_mr_auto():
    mtz = "test/beta_blip.mtz"
    pdb_list = {
        "test/beta.pdb": "1",
        "test/blip.pdb": "1"
    }
    assert mr_auto(mtz, pdb_list, None)
