"""trim_chothia module docstring.

This module provides tools for extracting Fv regions from numbered antibody PDB files.

Numbered pdb files can be downloaded from sabdab (http://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/search/?all=true)
however these contain remarks, heteroatoms, antigen residues and antibody residues that are not in the Fv region. 
This module allows the Fv atoms to be parsed out, and written to a separate file.

Example:
    $ python3 scripts/trim_chothia.py my_pdb.pdb C B my_pdb_trimmed.pdb  
"""
import argparse
import pathlib
import re
from typing import List

data_path = pathlib.Path(__file__).parents[1]/'data'/'example_pdbs'

def read_pdb(path: pathlib.PosixPath) -> List:
    return path.read_text().splitlines()

class Atom:
    """utility class for accessing atom attributes from raw string"""
    
    def __init__(self, raw: str) -> None:
        self.raw = raw
        
    def position(self) -> int:
        return int(re.match('\d+', self.raw[23:27].strip()).group())

    @property
    def chain_id(self) -> str:
        return self.raw[21:22]

    @chain_id.setter
    def chain_id(self, value) -> None:
        raw_list = list(self.raw)
        raw_list[21] = value
        self.raw = ''.join(raw_list)

def get_atoms(lines: List[str]) -> List[Atom]:
    return [Atom(line) for line in lines if 'ATOM' in line]

class Fv:
    """utility class for extracting Fv region from PDB and writing to disk
    
    Args:
        hid (str): the Fv heavy chain id
        lid (str): the Fv light chain id
    
    """
    heavy_region = list(range(1,114))
    light_region = list(range(1,108))

    def __init__(self, atoms):
        self.atoms = atoms

    def write_to_pdb(self, path):
        with open(path, 'w') as w:
            for atom in self.atoms:
                w.write(atom.raw + '\n')

    def reset_chain_id(self, old, new):
        # useful for changing heavy and light chain IDs to H and L respectively 
        # so they can be used as input to abangle
        for atom in self.atoms:
            if atom.chain_id == old:
                atom.chain_id = new
            else:
                continue

    @classmethod
    def from_atoms(cls, atoms, hid, lid):
        fv = []
        for atom in atoms:
            if atom.chain_id == hid and atom.position() in cls.heavy_region:
                fv.append(atom)
            elif atom.chain_id == lid and atom.position() in cls.light_region:
                fv.append(atom)
            else:
                continue
        return cls(fv)

parser = argparse.ArgumentParser(description='Process pdb files containing regions')
parser.add_argument('file', type=str, help='pdb file containing unparsed structure')
parser.add_argument('h', type=str, help='heavy chain id')
parser.add_argument('l', type=str, help='light chain id')
parser.add_argument('o', type=str, help='name of the output file')
args = parser.parse_args()

if __name__ == '__main__':

    pdb_path = data_path/args.file
    pdb = read_pdb(pdb_path)
    atoms = get_atoms(pdb)
    fv = Fv.from_atoms(atoms, args.h, args.l)
    for new_cid, old_cid in zip(['H', 'L'],[args.h, args.l]):
        fv.reset_chain_id(old_cid, new_cid)
    fv.write_to_pdb(data_path/args.o)


