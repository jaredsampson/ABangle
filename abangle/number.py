"""Module contains logic for locating and renumbering Fv chains in PDB files"""

import argparse
import pathlib
import re
from collections import defaultdict, namedtuple
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union, TextIO
from anarci import anarci
from dataclasses import dataclass

aa_codes_3to1 = {
        'ARG': 'R','HIS': 'H','LYS': 'K','ASP': 'D','GLU': 'E','SER': 'S', 'THR': 'T','ASN': 'N','GLN': 'Q','CYS': 'C',
        'GLY': 'G','PRO': 'P','ALA': 'A','VAL': 'V','ILE': 'I','LEU': 'L','MET': 'M','PHE': 'F','TYR': 'Y','TRP': 'W'
    }

@dataclass
class AtomRecord:
    """Class for accessing, setting and formatting atom records in a pdb file."""
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

    def reset_attribute(self, attr: str, replacement: Any) -> None:
        """Method allows you to edit the attributes of an atom corresponding to a single line in a pdb file.

        Args:
            attr (str): The attribute to be edited. 
            replacement (Any): The replacement value.
        """
        assert attr in self.__dict__.keys(), 'Attribute could not be found'
        correct_type = type(getattr(self, attr))
        if not isinstance(replacement, correct_type):
            try:
                replacement = correct_type(replacement)
            except ValueError:
                raise ValueError(f'Could not cast {attr} to {correct_type}')
        setattr(self, attr, replacement)
        
    @property
    def pdb_formatted_string(self) -> str:
        """Method for correctly formatting the attributes of an atom record for writing to a pdb file.

        Returns:
            str: A string representation of an atom that conforms to the pdb file format guide found
            here http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM. 
        """
        return "{:6s}{:5d}  {:<3s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n".format(*self.__dict__.values()) #Code borrowed from https://cupnet.net/pdb-format/

    @classmethod
    def from_string(cls, raw: str) -> None:
        """Method for constructing an AtomRecord object from a line in a pdb file.
        Can be used in list comp i.e. '[AtomRecord.from_string(line) for line in pdb_file_handle.readlines()].

        Args:
            raw (str): A raw string representing one line in a pdb file.

        Returns:
            AtomRecord: An AtomRecord object.
        """

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

class AtomRecordCollection:
    """Container for AtomRecord instances that allows convenient editing of multiple records simultaneously."""
    def __init__(self, name: str, records: List[AtomRecord]) -> None:
        self.name = name
        self.records = records

    def write_to_pdb(self, path: pathlib.Path) -> None:
        """Method for writing atom records to a pdb file with correct formatting.

        Args:
            path (pathlib.Path): Output file path.
        """

        with open(path, 'w') as f:
            for atom in self.records:
                f.write(atom.pdb_formatted_string)

    def _get_new_numbering(self) -> Tuple:
        """Method generates new numbering for antibody Fv chains based on the chothia numbering scheme. 
        Info on numbering can be found at https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6198058/.

        Raises:
            NotImplementedError: Raised if anarci failed to find H and L Fv chains.

        Returns:
            Tuple: (
                keys: Dictionary mapping concatenated chain, number, insertion code to residue (e.g. {'A_1_': 'T', ..., 'A_100_B': 'Y', ...}),
                numbering: Nested data structure containing number mapped to residue (e.g. [((1, ' '), 'T'), ..., ((100, 'B'), 'Y'), ...]),
                details: List of dictionaries containing details of the alignment e.g. start and end index of the chain in the query sequence, chain identifier etc
            )
        """
        keys, residues = list(zip(*self.indexed_sequence))
        sequence = ''.join(residues)
        numbering, details, _ = anarci([(self.name, sequence)], scheme = 'chothia')
        details, numbering = details[0], numbering[0] 
        if not any(numbering):
            raise NotImplementedError('Chains could not be renumbered')
        assert len(numbering) == 2, 'Only one chain could be renumbered'
        return keys, numbering, details

    def _map_old_to_new_numbering(self) -> Dict:
        """Method provides a dictionary to ensure new numbering is correctly applied to the appropriate AtomRecord. 
        A unique key is generated for each residue containing its chain, number and insertion code. 
        The start and end indices of the Fv chains are used to subset these keys and align them with the new numbering. 
        The resulting dictionary can be used to look up the correct new number for each atom.

        Raises:
            ValueError: Raised if chain is not a recognised antibody chain.

        Returns:
            Dict: Dictionary mapping concatenated chain, number, insertion code to residue (e.g. {'A_1_': 'T', ..., 'A_100_B': 'Y', ...}).
        """
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


    def renumber_residues(self) -> None:
        """Iterates through AtomRecord object and creates unique residue key. Key is then used to look up the.
        new chain id, number and insertion code in new numbering scheme. New numbering is applied and records.
        not in the Fv are dripped.
        """
        key_mapping = self._map_old_to_new_numbering()
        for atom in self.records:
            key = self._create_key(atom.chain_identifier, atom.residue_sequence_number, atom.insertion_code)
            
            if key in key_mapping:
                chain, number, insertion_code = key_mapping[key].split('_')
                replacement_attrs = zip(
                    ['chain_identifier', 'residue_sequence_number', 'insertion_code'], 
                    [chain, number, insertion_code]
                )
                
                for attr, replacement in replacement_attrs:
                    atom.reset_attribute(attr, replacement)
        
        self.records = [atom for atom in self.records if atom.chain_identifier in ('H', 'L')]

    @property
    def indexed_sequence(self) -> List[Tuple]:
        """Method iterates through alpha carbon records to create mapping between residue key: residue e.g. [..., ('A_100_A', 'G'), ('A_100_B', 'H')].

        Returns:
            List[Tuple]: List of residues mapped to their unique residue key.
        """

        return [
            (self._create_key(atom.chain_identifier, atom.residue_sequence_number, atom.insertion_code), aa_codes_3to1[atom.residue_name])
            for atom in self.records
            if atom.atom_name == 'CA' # only use alpha carbon
        ]

    @classmethod
    def from_file(cls, path: pathlib.Path) -> None:
        """Class method for reading in a collection of AtomRecords from a pdb file path and returning an AtomRecordCollection instance.

        Args:
            path (pathlib.Path): Path to pdb file e.g. 'data/pdb_files/1U8L.pdb'.

        Returns:
            AtomRecordCollection.
        """
        
        if not isinstance(path, pathlib.Path):
            path = pathlib.Path(path)

        return cls(
            path.stem, [
            AtomRecord.from_string(line) 
            for line in path.read_text().splitlines() 
            if line.startswith('ATOM')
        ])
    
    @classmethod
    def from_records(cls, name: str, records: List[AtomRecord]) -> None:
        """Class method for creating an AtomRecordCollection instance from a list of AtomRecords.

        Args:
            name (str): Name of the structure the AtomRecords belong to.
            records (List[AtomRecord]): A list of AtomRecord instances belonging to a single pdb structure.

        Returns:
            AtomRecordCollection.
        """
        return cls(name, records)

    @staticmethod
    def _create_key(chain: str, number: int, insertion_code: str) -> str:
        """Method for creating a key for each AtomRecord that can be used to look up the correct number when renumbering Fv chains.

        Args:
            chain (str): Chain identifier e.g. 'A'.
            number (int): Residue number e.g. 100.
            insertion_code (str): e.g. 'B' or ''.

        Returns:
            str: Unique residue key e.g. 'A_100_B'.
        """
        return '_'.join([chain, str(number), insertion_code])
    
    def __repr__(self) -> str:
        return f'{len(self.records)} records from structure {self.name}'

if __name__ == '__main__':
    pass