from abangle.calculate import find_angles
from abangle.number import *
import pytest
from Bio.PDB.Structure import Structure
from Bio.PDB.PDBExceptions import BiopythonWarning
import warnings

@pytest.fixture
def structure(shared_datadir):
    return get_structure(shared_datadir / '1u8l.pdb')

def test_get_structure(shared_datadir):
    assert isinstance(get_structure(shared_datadir / '1u8l.pdb'), Structure)

def test_get_structure_sequence(structure):
    assert True 

def test_find_angles_does_not_emit_biopython_warnings(shared_datadir):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        find_angles(shared_datadir / "1u8l.pdb")

    assert not any(issubclass(item.category, BiopythonWarning) for item in caught)
