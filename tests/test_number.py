from abangle.number import *
import pytest
import tempfile

@pytest.fixture
def pdb_string(datadir):
    return (datadir / 'pdb_str.txt').read_text()

def test_AtomRecord_constructor_from_string(pdb_string):
    atom_record = AtomRecord.from_string(pdb_string.splitlines()[0])
    assert atom_record.atom == 'ATOM'
    assert atom_record.atom_serial_number == 3338
    assert atom_record.atom_name == 'N'
    assert atom_record.alternate_location_indicator == ''
    assert atom_record.residue_name == 'SER'
    assert atom_record.chain_identifier == 'B'
    assert atom_record.residue_sequence_number == 215
    assert atom_record.insertion_code == ''
    assert atom_record.x_coordinate == 6.364
    assert atom_record.y_coordinate == 9.072
    assert atom_record.z_coordinate == 33.809
    assert atom_record.occupancy == 1.00
    assert atom_record.element_symbol == 'N'  

def test_AtomRecord_reset_attribute_chain_id(pdb_string):
    atom_record = AtomRecord.from_string(pdb_string.splitlines()[0])
    assert atom_record.chain_identifier == 'B'
    atom_record.reset_attribute('chain_identifier', 'H')
    assert atom_record.chain_identifier == 'H'

def test_AtomRecord_reset_attribute_residue_sequence_number(pdb_string):
    atom_record = AtomRecord.from_string(pdb_string.splitlines()[0])
    assert atom_record.residue_sequence_number == 215
    atom_record.reset_attribute('residue_sequence_number', 216)
    assert atom_record.residue_sequence_number == 216

def test_AtomRecord_reset_attribute_wrong_type(pdb_string):
    atom_record = AtomRecord.from_string(pdb_string.splitlines()[0])
    with pytest.raises(ValueError):
        assert atom_record.residue_sequence_number == 215
        atom_record.reset_attribute('residue_sequence_number', 'B')

def test_AtomRecord_pdb_formatted_string(pdb_string):
    atom_lines = [line for line in pdb_string.splitlines() if line.startswith('ATOM')]
    records = [AtomRecord.from_string(line) for line in atom_lines]
    assert all([
        record.pdb_formatted_string.strip() == line.strip() # watch out for whitespace/\n chars at the end of string
        for record, line in zip(records, atom_lines)
    ])

def test_AtomRecordCollection_from_records(pdb_string):
    records = [
        AtomRecord.from_string(line) 
        for line in pdb_string.splitlines() 
        if line.startswith('ATOM')
    ]
    record_collection = AtomRecordCollection.from_records('1PDB', records)
    assert len(record_collection.records) == 21
    assert record_collection.name == '1PDB'

def test_AtomRecordCollection_indexed_sequence(datadir):
    atom_record_collection = AtomRecordCollection.from_file(datadir / 'pdb_str.txt')
    assert atom_record_collection.indexed_sequence == [('B_215_', 'S'), ('B_216_', 'C'), ('C_1_', 'D')]

def test_AtomRecordCollection_create_key():
    assert AtomRecordCollection._create_key('A', 105, 'A') == 'A_105_A'
    assert AtomRecordCollection._create_key('B', 32, '') == 'B_32_'
    assert AtomRecordCollection._create_key('A', 1, '') == 'A_1_'
    assert AtomRecordCollection._create_key('L', 250, 'C') == 'L_250_C'

def test_AtomRecordCollection_write_to_pdb(tmp_path, pdb_string, datadir):
    atom_record_collection = AtomRecordCollection.from_file(datadir / 'pdb_str.txt')
    d = tmp_path / "sub"
    d.mkdir()
    p = d / "test.pdb"
    atom_record_collection.write_to_pdb(p)
    assert p.read_text().splitlines() == [
        line for line in pdb_string.splitlines() 
        if line.startswith('ATOM')
    ]

