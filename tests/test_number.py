from abangle.number import *
import pytest
from Bio.PDB.Structure import Structure

@pytest.fixture
def structure(shared_datadir):
    return get_structure(shared_datadir, '1U8L.pdb')

def test_get_structure(shared_datadir):
    assert isinstance(get_structure(shared_datadir, '1U8L.pdb'), Structure)

def test_get_structure_sequence(structure):
    assert True 