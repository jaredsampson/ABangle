import argparse
import pathlib
import re
from collections import defaultdict, namedtuple
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union, TextIO
from anarci import anarci
from dataclasses import dataclass

AlignmentDetails = list[list[dict[str, Any]]]

aa_codes_3to1 = {
        'ARG': 'R','HIS': 'H','LYS': 'K','ASP': 'D','GLU': 'E','SER': 'S', 'THR': 'T','ASN': 'N','GLN': 'Q','CYS': 'C',
        'GLY': 'G','PRO': 'P','ALA': 'A','VAL': 'V','ILE': 'I','LEU': 'L','MET': 'M','PHE': 'F','TYR': 'Y','TRP': 'W'
    }

@dataclass
class Atom:
    atom: str
    atom_serial_number: str
    atom_name: str
    alternate_location_indicator: str
    residue_name: str
    chain_identifier: str
    residue_sequence_number: str
    insertion_code: str
    x_coordinate: str
    y_coordinate: str
    z_coordinate: str
    occupancy: str
    temperature_factor: str
    element_symbol: str
    charge: str

    def reset_attribute(self, attr, replacement):
        assert attr in self.__dict__.keys(), 'Attribute could not be found'
        correct_type = type(getattr(self, attr))
        if not isinstance(replacement, correct_type):
            try:
                replacement = correct_type(replacement)
            except ValueError:
                print(f'Could not cast {attr} to {correct_type}')
        setattr(self, attr, replacement)
        
    @property
    def pdb_formatted_string(self):
        return "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n".format(*self.__dict__.values())

    @classmethod
    def from_string(cls, raw: str):
        assert raw.startswith('ATOM'), 'Line does not describe an atom'
        
        indices = {
            'atom': (0,6),
            'atom_serial_number': (6,11),
            'atom_name': (12,16),
            'alternate_location_indicator': (16,17),
            'residue_name': (17,20),
            'chain_identifier': (21,22),
            'residue_sequence_number': (22,26),
            'insertion_code': (26,27),
            'x_coordinate': (30,38),
            'y_coordinate': (38,46),
            'z_coordinate': (46,54),
            'occupancy': (54,60),
            'temperature_factor': (60,66),
            'element_symbol': (76,78),
            'charge': (78,80)
        }

        attrib_types = (
            str, int, str, str, str, str, int, str, float, float, float, float, float, str, str
        )
        
        return cls(**{
            key: new_type(raw[slice(*ind)].strip()) 
            for (key, ind), new_type in zip(indices.items(), attrib_types)
        })

class AtomList:
    def __init__(self, name, atoms: List[Atom]):
        self.name = name
        self.atoms = atoms

    def write_to_pdb(self, path: pathlib.Path):
        with open(path, 'w') as f:
            for atom in self.atoms:
                f.write(atom.pdb_formatted_string)

    def _get_new_numbering(self) -> Tuple:
        keys, residues = list(zip(*self.indexed_sequence))
        sequence = ''.join(residues)
        numbering, details, _ = anarci([(self.name, sequence)], scheme = 'chothia')
        details, numbering = details[0], numbering[0] 
        if not any(numbering):
            raise NotImplementedError('Chains could not be renumbered')
        assert len(numbering) == 2, 'Only one chain could be renumbered'
        return keys, numbering, details

    def _map_old_to_new_numbering(self):
        old_keys, numbering, details = self._get_new_numbering()
        key_map = []
        for num, det in zip(numbering, details):
            num = num[0]
            
            if det['chain_type'] in ('K', 'L', 'H'):
                chain = 'L' if det['chain_type'] in ('K', 'L') else 'H' 
            else:
                raise ValueError('Chain type not recognized')
            
            new_keys = [
                self._create_key(chain, str(n[0][0]), n[0][1]).strip()
                for n in num 
                if not n[1] == '-'
            ]
            key_map.extend(list(zip(old_keys[det['query_start']: det['query_end']], new_keys)))
        
        return dict(key_map)

    def _create_key(self, chain, number, insertion_code):
        return ('_'.join([chain, str(number), insertion_code]))

    def renumber(self):
        key_mapping = self._map_old_to_new_numbering()
        for atom in self.atoms:
            key = self._create_key(atom.chain_identifier, atom.residue_sequence_number, atom.insertion_code)
            
            if key in key_mapping:
                chain, number, insertion_code = key_mapping[key].split('_')
                replacement_attrs = zip(
                    ['chain_identifier', 'residue_sequence_number', 'insertion_code'], 
                    [chain, number, insertion_code]
                )
                
                for attr, replacement in replacement_attrs:
                    atom.reset_attribute(attr, replacement)
        
        self.atoms = [atom for atom in self.atoms if atom.chain_identifier in ('H', 'L')]

    @property
    def indexed_sequence(self) -> List[Tuple]:
        return [
            (self._create_key(atom.chain_identifier, atom.residue_sequence_number, atom.insertion_code), aa_codes_3to1[atom.residue_name])
            for atom in self.atoms
            if atom.atom_name == 'CA' # only use alpha carbon
        ]

    @classmethod
    def from_file(cls, path: pathlib.Path):
        return cls(
            path.stem, [
            Atom.from_string(line) 
            for line in path.read_text().splitlines() 
            if line.startswith('ATOM')
        ])
    
    @classmethod
    def from_atoms(cls, name, atoms: List[Atom]): #-> AtomList how to provide type hint for class in its own definition?
        return cls(name, atoms)

    def __repr__(self) -> str:
        return f'{len(self.atoms)} atoms from structure {self.name}'

if __name__ == '__main__':
    data_path = pathlib.Path().parent/'tests'/'data'
    pdb_original = (data_path/'1U8L.pdb').read_text().splitlines()
    atoms = [Atom.from_string(line) for line in pdb_original if line.startswith('ATOM')]
    atomlist = AtomList.from_file(data_path/'1U8L.pdb')
    atomlist.renumber()
    atomlist.write_to_pdb(data_path/'test.pdb')