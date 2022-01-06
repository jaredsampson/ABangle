import pathlib
import argparse
from typing import *
from anarci import anarci
from fastcore.basics import Str
from fastcore.xtras import dict2obj

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO, Select
from Bio.PDB.Structure import Structure
from Bio.PDB.Chain import Chain
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def get_structure(path: pathlib.Path) -> Structure:
    """Parses PDB file, returns a structure object.

    Args:
        data_path (pathlib.Path): Path to data
    """
    parser = PDBParser(PERMISSIVE=True, QUIET=True)
    return parser.get_structure(path.stem, path)

def get_structure_sequences(structure: Structure) -> Dict[str, SeqRecord]:
    """Gets the sequences of each chain as a list of SeqRecord objects 
    """
    return {
        seq.id[-1]: seq # 1U8L:A[-1] -> A
        for seq in SeqIO.PdbIO.AtomIterator(
            structure.id, structure=structure)
        }

def parse_chain_type(details: Dict[str, Any]) -> str:
    """Takes a one letter chain identifier and converts it to either H, L or raises error

    Args:
        details (Dict): A dictionary containing details of the antibody alignment, chain
        type will be either H (heavy), L (lambda) or K (kappa). Abangle expects H (heavy)
        and L (light) chain identifiers, so function will convert L or K -> L     

    Returns:
        str: A chain identifier 
    """
    chain_type = details['chain_type']
    if chain_type in ('L', 'K'):
        return 'L'
    elif chain_type == 'H':
        return 'H'
    else:
        raise ValueError(f'chain type {chain_type}')
    
def number_sequences(sequences: Dict[str, SeqRecord], scheme: str = 'chothia') -> Dict[str, Dict[str, Any]]:
    """Takes a list of chain: sequence mappings and aligns the sequences to get their chain identifer and numbering
    The numbering is then returned as a dictionary mapping query chain name to Fv chain, span and

    Args:
        sequences (Dict[str, SeqRecord]): 
        scheme (str, optional): Numbering scheme to use. Can be one of: imgt, chothia, kabat or martin. Defaults to 'chothia'. 

    Returns:
        Dict[str, Dict[str, Any]]: e.g. {'A': {'Chain': 'H', 'span': slice(0, 107, None), numbering: (' ', 1', ' ')}}
    """
    input = [
        (chain, str(seq.seq)) 
        for chain, seq in sequences.items()
    ]
    
    numbering, details, _ = anarci(input, scheme = scheme)
    
    numbering = [
        [(' ', res_id[0][0], res_id[0][1]) for res_id in num[0][0]
        if res_id[1] != '-'] # ungap the sequence
        for num in numbering 
        if num
    ]
    
    details = dict2obj([det[0] for num, det in zip(numbering, details) if num]) # details wrapped in list

    return dict2obj({
        det.query_name: 
        {'chain': parse_chain_type(det),'span': slice(det.query_start, det.query_end), 'numbering': num}
        for num, det in zip(numbering, details)
    })

def renumber_structure(structure: Structure, numbering: Dict) -> None:
    """Takes a numbering dictionary and a structure and updates the residue IDs with the new numbering 

    Args:
        structure (Structure): 
        numbering (Dict):
    """
    for name, numbering in numbering.items():
        res_ids = [res.id for res in structure[0][name].get_residues()][numbering.span]
        
        for old, new in zip(res_ids, numbering.numbering):
            structure[0][name][old].id = new
        
        structure[0][name].id = numbering.chain

class FvSelect(Select):
    """Sublassess Select so that only residues in the Fv will be written to disk

    Args:
        Select ([type]): [description]
    """
    def accept_residue(self, residue):
        
        hetcode, seqid, icode = residue.id
        chain = residue.parent.id
        
        if chain == 'H' and seqid in range(113):
            return True
        elif chain == 'L' and seqid in range(107):
            return True
        else:
            return False

def write_pdb(path: str, structure: Structure) -> None:
    """Writes the Fv portion of a pdb to disk
    """
    io = PDBIO()
    io.set_structure(structure)
    io.save(path, select=FvSelect())

def main(argv=None):
    parser = argparse.ArgumentParser(description='Renumber the Fv portion of a pdb file')
    parser.add_argument('infile')
    parser.add_argument('outfile')
    args = parser.parse_args(argv)
    structure = get_structure(pathlib.Path(args.infile))
    sequences = get_structure_sequences(structure)
    numbering = number_sequences(sequences)
    renumber_structure(structure, numbering)
    write_pdb(args.outfile, structure)

if __name__ == '__main__':
    exit(main())
    