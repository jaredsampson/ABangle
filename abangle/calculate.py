"""
DESCRIPTION	
	
	abangle.calculate.py
	A module to calculate the orientation between the VH and VL domains in an Antibody Fv regions.

OUTPUT

	~ Angles of the Fv regions in the files inputted from command line interface.

AUTHOR
	
	2012
	James Dunbar - Oxford protein informatics group.
	james.dunbar@dtc.ox.ac.uk
	www.stats.ox.ac.uk/~dunbar/
		
	Supervisors
	Prof C.Deane - Oxford protein informatics group.
	Dr Angelika Fuchs (Roche) and Dr Jiye Shi (UCB Celltech)\n
"""
import subprocess
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.Structure import Structure
import numpy
import math
import sys
import os
import tempfile
import urllib.request, urllib.parse, urllib.error
import itertools
import pathlib
from Bio.PDB.PDBParser import PDBParser
from abangle import dataIO
from abangle import number as num
from abangle.constants import coreset, centroids, aa_code_dict
from typing import List
from dataclasses import dataclass

path = pathlib.Path(__file__).parent
data_path = path.parent/'data'

#########################
# Calculation functions #
#########################

def create_coresets(path):
    # create paths
    name = path.stem[:4]
    Houtpath = str((path.parent/f'{path.stem}_Hcoreset').with_suffix('.pdb'))
    Loutpath = str((path.parent/f'{path.stem}_Lcoreset').with_suffix('.pdb'))
    
    # read in structure
    parser = PDBParser()
    structure = parser.get_structure(name, path)
    
    # save H and L coresets to separate files
    io = PDBIO()
    io.set_structure(structure)
    io.save(Houtpath, select = num.HCoresetSelect())
    io.save(Loutpath, select = num.LCoresetSelect())

    # return the filehandles
    return Houtpath, Loutpath

def mapvectors(fname, PAPS_def=False):
    """Maps the reference frames (planes) onto to VH and VL domains of an Fv structure (fname is chothia numbered pdb file
    with VH as H chain and VL as L chain. PAPS_def means use the same definition of C that Abhinandan  and Martin did when calculating
    their torsion angle (makes HL should be the same as their packing angle as defined in authors' paper)"""
    # Get transformation matrices by aligning the core of the domains
    Hf, Lf = create_coresets(fname)
    
    uL = align(os.path.join(data_path, "consensus_L.pdb"), Lf)
    uH = align(os.path.join(data_path, "consensus_H.pdb"), Hf)
    os.remove(Hf)
    os.remove(Lf)

    if PAPS_def:
        # The centroids of interface residues.
        cH = centroids['heavy'][2]
        cL = centroids['light'][2]
    else:
        # The minimally varying centroid vector is at. As calculated.:
        cH = [
            -10 * 0.5 * centroids['heavy'][0][i] + 1 * 0.5 * centroids['heavy'][1][i] + centroids['heavy'][2][i] for i in range(3)
        ]
        cL = [
            6 * 0.5 * centroids['light'][0][i] - 2 * 0.5 * centroids['light'][1][i] + centroids['light'][2][i] for i in range(3)
        ]

    # Define the plane vectors from the centroid point
    # On VL domain
    L1 = [cL[i] + centroids['light'][0][i] for i in range(3)]
    L2 = [cL[i] + centroids['light'][1][i] for i in range(3)]

    # On VH domain
    H1 = [cH[i] + centroids['heavy'][0][i] for i in range(3)]
    H2 = [cH[i] + centroids['heavy'][1][i] for i in range(3)]

    # Do the transfomation onto the
    Lpoints = list([transform(x, uL) for x in (cL, L1, L2)])
    Hpoints = list([transform(x, uH) for x in (cH, H1, H2)])

    return Lpoints, Hpoints

def align(file1, file2):
    """Aligns file1 to file2 using tmalign and returns the transformation matrix"""
    # Temp file for the matrix for latest versions of TMalign
    mtmpfd, mtmp = tempfile.mkstemp(".txt", "matrix")
    os.close(mtmpfd)
    # Align file1 to file2 using TMalign

    try:
        subpr = subprocess.Popen(
            ["TMalign", file1, file2, "-m", mtmp],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        TMresult = subpr.communicate()
    except OSError:
        # If is not found, point to webpage for installation.
        raise Exception(
            "Cannot execute TMalign. Please install and ensure it is in your path.\nTMalign can be downloaded from:\n"
            "http://zhanglab.ccmb.med.umich.edu/TM-align/\n"
            "Reference: Y. Zhang and J. Skolnick, Nucl. Acids Res. 2005 33, 2302-9\n"
        )

    # Parse the output of TMalign. Some versions don't output the matrix. -m option is needed. Does not affect versions which don't need it.
    
    result = TMresult[0].decode('utf-8').split("\n")
    attempt = 0
    while 1:
        try:
            i = 0
            while 1:
                if result[i].upper().startswith(" -------- ROTATION MATRIX"):
                    # Grab transformation matrix
                    u = []
                    u.append(list(map(float, result[i + 2].split()[1:])))
                    u.append(list(map(float, result[i + 3].split()[1:])))
                    u.append(list(map(float, result[i + 4].split()[1:])))
                    break
                else:
                    i += 1
            break
        except IndexError:
            try:
                if not attempt:
                    ftmp = open(mtmp)
                    result = ftmp.readlines()
                    ftmp.close()
                    attempt = 1
            except IOError:
                break

    if os.path.exists(mtmp):
        os.remove(mtmp)

    # Return the transformation matrix
    try:
        return u
    except NameError:
        raise Exception(
            "TMalign alignment file not in an expected format, check output gives rotation matrix (or with -m option )\n"
        )

def transform(coords, u):
    """Transforms coords by a matrix u. u is found using tmalign"""
    # Ensure coordinates are of type float
    coords = list(map(float, coords))
    # Do transformation
    X = u[0][0] + u[0][1] * coords[0] + u[0][2] * coords[1] + u[0][3] * coords[2]
    Y = u[1][0] + u[1][1] * coords[0] + u[1][2] * coords[1] + u[1][3] * coords[2]
    Z = u[2][0] + u[2][1] * coords[0] + u[2][2] * coords[1] + u[2][3] * coords[2]

    # Return transformed coordinates
    return [X, Y, Z]

def normalise(vec):
    mag = (sum(list([x ** 2 for x in vec]))) ** 0.5
    return list([x / mag for x in vec])

def angles(fname):
    """Calculate the orientation measures for the structure in fname"""
    # Map the vectors on the Heavy and Light domains of the structure
    
    Lpoints, Hpoints = mapvectors(fname)

    # Create vectors with which to calculate angles between.
    C = normalise([Hpoints[0][i] - Lpoints[0][i] for i in range(3)])
    Cminus = list([-1 * x for x in C])
    L1 = normalise([Lpoints[1][i] - Lpoints[0][i] for i in range(3)])
    L2 = normalise([Lpoints[2][i] - Lpoints[0][i] for i in range(3)])
    H1 = normalise([Hpoints[1][i] - Hpoints[0][i] for i in range(3)])
    H2 = normalise([Hpoints[2][i] - Hpoints[0][i] for i in range(3)])
    dc = (
        sum([x ** 2 for x in [Hpoints[0][i] - Lpoints[0][i] for i in range(3)]])
        ** 0.5
    )

    # Projection of the L1 and H1 vectors onto the plane perpendicular to the centroid vector.
    n_x = numpy.cross(L1, C)
    n_y = numpy.cross(C, n_x)

    tmpL_ = normalise([0, numpy.dot(L1, n_x), numpy.dot(L1, n_y)])
    tmpH_ = normalise([0, numpy.dot(H1, n_x), numpy.dot(H1, n_y)])

    # HL is the angle between the L1 and H1 vectors looking down the C vector (the centroid vector)
    HL = math.acos(numpy.dot(tmpL_, tmpH_))
    HL = HL * (180.0 / math.pi)

    # Find direction by computing cross products
    if numpy.dot(numpy.cross(tmpL_, tmpH_), [1, 0, 0]) < 0:
        HL = -HL

    # LC1 angle is the angle between the L1 and C vectors
    LC1 = math.acos(numpy.dot(L1, C))
    LC1 = LC1 * (180.0 / math.pi)

    # HC1 angle is the angle between the H1 and C vectors
    HC1 = math.acos(numpy.dot(H1, Cminus))
    HC1 = HC1 * (180.0 / math.pi)

    # LC2 angle is the angle between the L2 and C vectors
    LC2 = math.acos(numpy.dot(L2, C))
    LC2 = LC2 * (180.0 / math.pi)

    # HC2 angle is the angle between the H2 and C vectors
    HC2 = math.acos(numpy.dot(H2, Cminus))
    HC2 = HC2 * (180.0 / math.pi)

    # Return the angles and the separation distance.
    return dict(
        list(zip(["HL", "HC1", "LC1", "HC2", "LC2", "dc"], [HL, HC1, LC1, HC2, LC2, dc]))
    )

def get_loop_length(seq, loop):
    
    if loop == 'L1':
        residues = ["L24", "L25", "L26", "L27", "L28", "L29", "L30", "L31", "L32", "L33", "L34"]
    elif loop == 'H3':
        residues = ["H95", "H96", "H97", "H98", "H99", "H100", "H101", "H102"]
    else:
        raise ValueError('Loop not recognized')
    
    n = 0
    for res, aa in seq.items():
        if aa != "-":
            if res[-1].isdigit():
                if res in residues:
                    n += 1
            elif res[-1].isalpha():
                if res[:-1] in residues:
                    n += 1
            else:
                continue
    
    return n

def validate_angles(file, HC2, HC1, LC2, LC1, dc, HL, rel_tol=1e-2):
    name = file[:4]
    angle_keys = ['HC2', 'HC1', 'LC2', 'LC1', 'dc', 'HL']
    true_angles = [HC2, HC1, LC2, LC1, dc, HL]
    new_angles = angles(examples/file)
    assert all([math.isclose(new_angles[key], angle, rel_tol=rel_tol) for key, angle in zip(angle_keys, true_angles)])

if __name__ == '__main__':
    
    examples = data_path.parent/'data'/'example_pdbs'

    angles_4KQ3 = {'HC2':114.97, 'HC1':71.58, 'LC2':83.15, 'LC1':119.49, 'dc':16.00, 'HL':-61.10}
    validate_angles('4kq3_chothia_Fv.pdb', **angles_4KQ3)

    angles_2ATK = {'HL': -54.09, 'HC1': 70.76, 'HC2': 114.09, 'LC1': 123.29, 'LC2': 83.42, 'dc': 16.48}
    validate_angles('2atk_chothia_Fv.pdb', **angles_2ATK)

    angles_1U8L = {'HL': -56.41, 'HC1': 68.73, 'HC2': 118.87, 'LC1': 123.11, 'LC2': 82.99, 'dc': 15.88}
    validate_angles('1u8l_chothia_Fv.pdb', **angles_1U8L)
