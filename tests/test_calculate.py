from typing import List
from abangle.calculate import *
import argparse
import pytest
from pathlib import Path
from dataclasses import dataclass

path = Path(__file__).parent
data_path = path.parent/'data'

@dataclass
class PDBargs:
    scfv: List = None
    q: bool = True
    usernumbered: bool = True

@dataclass
class GetAnglesargs:
    i: str 
    store: bool = False
    q: bool = True

@pytest.fixture
def pdb():
    args = PDBargs()
    name = '4KQ3'
    filepath = data_path/'example_pdbs'/'4KQ3_abnum.pdb'
    return PDB(name, filepath, args)

def test_parse_pdb(pdb):
    pdb.parse_pdb()
    assert list(pdb.Chains.keys()) == [('', 'H'), ('', 'L')]
    assert pdb.ABchains == [('', 'H'), ('', 'L')]

def test_pair_chains(pdb):
    pdb.parse_pdb()
    pdb.pair_chains()
    fvs = pdb.Fvs[0]
    seq = fvs.sequence
    assert fvs.name == '4KQ3_HL'
    assert fvs.c1.FVATOMS[0] == 'ATOM      1  N   VAL H   2      28.598  -0.497   9.629  1.00 31.59           N  '
    assert fvs.c2.FVATOMS[0] == 'ATOM   1804  N   ASP L   1      24.272  29.366   0.981  1.00 46.85           N  '
    assert all([seq['H2'] == 'V', seq['H60'] == 'A', seq['L93'] == 'D', seq['H86'] =='D', seq['H100'] == 'F'])

def test_GetAngles():
    path = data_path/'example_pdbs'/'4KQ3_abnum.pdb'
    args = GetAnglesargs(path)
    return GetAngles(args)





    



