from abangle.number import *
import pytest
from Bio.PDB.Structure import Structure

@pytest.fixture
def structure(shared_datadir):
    return get_structure(shared_datadir / '1u8l.pdb')

def test_get_structure(shared_datadir):
    assert isinstance(get_structure(shared_datadir / '1u8l.pdb'), Structure)

def test_get_structure_sequence(structure):
    assert True 
