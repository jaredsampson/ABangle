"""
DESCRIPTION	
	
	abangle.calculate.py
	A module to calculate the orientation between the VH and VL domains in 
        an Antibody Fv regions.

OUTPUT

	~ Angles of the Fv regions in the files inputted from command 
        line interface.

AUTHOR
	
	2012
	James Dunbar - Oxford protein informatics group.
	james.dunbar@dtc.ox.ac.uk
	www.stats.ox.ac.uk/~dunbar/
		
	Supervisors
	Prof C.Deane - Oxford protein informatics group.
	Dr Angelika Fuchs (Roche) and Dr Jiye Shi (UCB Celltech)
"""
from Bio.PDB.Superimposer import Superimposer
from Bio.PDB.PDBParser import PDBParser

import json
import numpy as np
import math
import pathlib
from typing import List
from collections import namedtuple

from abangle.number import (
        get_structure_sequences,
        number_sequences,
        renumber_structure
        )

path = pathlib.Path(__file__).parent
data_path = path.parent/'data'

# load in coreset residue number dictionary
with open(data_path/'coresets.json') as f:
    coresets = json.load(f)

# Principal components computed on heavy chain consensus structure
pcH = np.loadtxt(data_path/'principal_components_heavy.csv')
pcL = np.loadtxt(data_path/'principal_components_light.csv')

def get_coreset_atoms(structure, chain, coresets):
    """Retrieves a list of Atom objects corresponding to the CA carbon of
    the coreset residues"""
    
    coreset_atoms = [
        atom
        for atom in structure.get_atoms()
        if atom.get_name() == 'CA'
        if atom.parent.id[1] in coresets[chain]
        if atom.parent.parent.id == chain
    ]
    
    return coreset_atoms

def align(coresets, structure, chain, consensus):
    """Computes rotation and translation matrices that can be used to map 
    vectors calculated on consensus structure onto query structure."""

    coreset_atoms = get_coreset_atoms(structure, chain, coresets)
    consensus_atoms = list(consensus.get_atoms())

    si = Superimposer()
    si.set_atoms(coreset_atoms, consensus_atoms)
    rot, tran = si.rotran
    return rot, tran

def transform(vector, rot, tran): 
    """Performs transformation of vector position given rotation and 
    transformation matrices."""
    return vector @ rot + tran

# Object holds the vectors that sit on a single plame
Points = namedtuple('Points', ['C', 'V1', 'V2']) 

def map_vectors(fname, chain, pcs, PAPS_def=False):
    """Maps the reference frames (planes) onto to VH and VL domains of an Fv structure (fname is chothia numbered pdb file
    with VH as H chain and VL as L chain.
    """
    # coefs to map centroids onto plane formed by PC1 and PC2 
    coefs = np.array([-5, 0.5, 1]) if chain == 'H' else np.array([3, -1, 1])
    
    # calculate the minimally varying centroid vector
    C = coefs @ pcs

    # Define the plane vectors from the centroid point
    V1 = C + pcs[0]
    V2 = C + pcs[1]

    parser = PDBParser(PERMISSIVE=True, QUIET=True)
    structure = parser.get_structure(fname.stem, fname)
    consensus = parser.get_structure(chain, data_path/f'consensus_{chain}.pdb')

    sequence = get_structure_sequences(structure)
    numbering = number_sequences(sequence, scheme = 'chothia')
    renumber_structure(structure, numbering)

    # Get transformation matrices by aligning the core of the domains
    rot, tran = align(coresets, structure, chain, consensus)

    # Do the transfomation onto the
    points = [transform(v, rot, tran) for v in (C, V1, V2)]

    return Points(*points)

def as_unit_vector(vec): 
    """Divides a vector by its length to give a vector of length 1 
    (unit vector)"""
    return vec / np.linalg.norm(vec)

def compute_angle(vec1, vec2): 
    """Computes angle between two vectors"""
    return np.arccos(vec1 @ vec2)

def find_angles(fname):
    """Calculate the orientation measures for the structure in fname"""
    
    RADIAN = 180.0 / math.pi
    
    # HL is the angle between the L1 and H1 vectors looking down the C vector (the centroid vector)
    # Map the vectors on the Heavy and Light domains of the structure
    Lpoints, Hpoints = [
        map_vectors(fname, chain, pcs) 
        for chain, pcs in zip(['L', 'H'], [pcL, pcH])
    ]

    C = as_unit_vector(Hpoints.C - Lpoints.C)
    L1 = as_unit_vector(Lpoints.V1 - Lpoints.C)
    L2 = as_unit_vector(Lpoints.V2 - Lpoints.C)
    H1 = as_unit_vector(Hpoints.V1 - Hpoints.C)
    H2 = as_unit_vector(Hpoints.V2 - Hpoints.C)

    dc = np.linalg.norm(Hpoints.C - Lpoints.C)

    # Projection of the L1 and H1 vectors onto the plane perpendicular to the centroid vector.
    n_x = np.cross(L1, C)
    n_y = np.cross(C, n_x)

    tmp_L = as_unit_vector([0, L1 @ n_x, L1 @ n_y])
    tmp_H = as_unit_vector([0, H1 @ n_x, H1 @ n_y])

    HL = compute_angle(tmp_L, tmp_H) * RADIAN

    # Find direction by computing cross products
    if np.cross(tmp_L, tmp_H) @ [1, 0, 0] < 0:
        HL = -HL
 
    LC1, HC1, LC2, HC2 = [
        compute_angle(vec1, vec2) * RADIAN 
        for vec1, vec2 in zip(
            [L1, H1, L2, H2], 
            [C, -C, C, -C]
        )
    ]

    # Return the angles and the separation distance.
    return {
        name: angle
        for name, angle in zip(
            ["HL", "HC1", "LC1", "HC2", "LC2", "dc"], 
            [HL, HC1, LC1, HC2, LC2, dc]
        )
    }

