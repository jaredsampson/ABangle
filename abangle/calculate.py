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

try:
    savedpath = open(os.path.join(path, "config/userdatapath.txt")).readline().strip()
    if savedpath and not os.path.exists(savedpath):
        raise IOError
    else:
        user_datapath = savedpath
except IOError:
    savedpath = ""
    user_datapath = ""

# Get the path to the antibody numbering program. (It may also be 'online' if the user has set this as default in the set up)
try:
    AnnotationProgPath = (
        open(os.path.join(path, "config/AnnotationProgPath.txt")).readline().strip()
    )
except IOError:
    AnnotationProgPath = ""

# Read in known sequences
Sequences, Residues = dataIO.load(
    os.path.join(data_path, "Sequences.dat"), header=True, rownames=0, conv=False
)

class PDB:
    """A class to describe PDB files. Parses the file and locates Fv regions.
    Should handle any structure with VH and VL regions in different chains. Use the scfv option to tell abnum to look for VH and VL regions on the same chain.
    Multiple models in the same file will also be handled. Note that abangle was developed using X-ray xtal structures only"""

    distance_cutoff = 20 # angstrom threshold for pairing heavy and light chains

    def __init__(self, name, filepath, args, verbose = False, usernumbered = True):
        self.name = name
        self.filepath = filepath
        self.args = args
        self.parser = PDBParser(PERMISSIVE=1)
        self.verbose = verbose
        self.usernumbered = usernumbered
        self.parse_pdb()  # Get Chain objects
        self.pair_chains()  # Pair the VH and VL chains using distance constraint

    def read_pdb(self):
        pdb_path = pathlib.Path(self.filepath)
        assert pdb_path.exists(), 'Could not find file, check the path is correct'
        return pdb_path.read_text().splitlines()
        
    def parse_pdb(self):
        # pdb_fo = io.StringIO(self.read_pdb())
        # pdb_lines = self.read_pdb()
        with open(self.filepath) as f:
            pdb_lines = f.readlines()
        
        model = ''
        chainlines = {}
        for line in pdb_lines:
            if line.startswith("MODEL"):
				# If we find a model id then use it for the following lines
                model = line.split()[1]
            elif line.startswith("ATOM"): # or line.startswith("HETATM"):
                c = line[21]
                try:
                    chainlines[ (model, c ) ].append( line )
                except:
                    chainlines[ (model, c ) ] = [ line ]

        # Create Chain objects for each chain in the file
        if self.verbose:
            if self.usernumbered:
                # Parse the PDB to get the correct annotation.
                sys.stdout.write(
                    f"Using the annotation of chain H and chain L in file {self.filepath}...\n"
                )
            elif self.args.online or AnnotationProgPath == "Use online":
                sys.stdout.write("Annotating using online server...\n")
            elif AnnotationProgPath:
                sys.stdout.write(
                    f"Annotating using {os.path.split(AnnotationProgPath)[1]} ...\n" 
                )

        self.Chains = {
            chain: Chain(chainlines[chain], chain, self.args) for chain in chainlines
        }
        self.ABchains = [
            chain for chain in self.Chains if self.Chains[chain].type in ["H", "L"]
        ]

        # TODO test with pdb in scfv format (both H and L are in single chain) with http://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/structureviewer/?pdb=6lfm
        # If the user has provided a chain id for a single chain fv then create a new chain object with the unannotated atoms
        if self.args.scfv and self.ABchains:
            if self.verbose:
                sys.stdout.write("Attempting split of scfv into VH and VL...\n")
            self.ScFv = {}
            for chain in self.Chains:
                if chain[1] in self.args.scfv:
                    self.ScFv[(chain[0], chain[1], 1)] = self.Chains[chain]
                    if self.verbose:
                        sys.stdout.write(
                            f"""Found model {self.Chains[chain].model} chain {self.Chains[chain].ident} V{self.Chains[chain].type} region.\n"""
                        )
                    # Annotate the remainder of the structure
                    self.ScFv[(chain[0], chain[1], 2)] = Chain(
                        self.Chains[chain].unannotated, chain, self.args
                    )
                    if self.ScFv[(chain[0], chain[1], 2)].type == "Antigen":
                        if self.verbose:
                            sys.stderr.write(
                                "Warning: Abnum failed to find V%s of your scFv %s chain %s\n"
                                % (
                                    ["H", "L"][
                                        int(
                                            not ["H", "L"].index(
                                                self.Chains[chain].type
                                            )
                                        )
                                    ],
                                    chain[0],
                                    chain[1],
                                )
                            )
                        del self.ScFv[(chain[0], chain[1], 1)]
                        del self.ScFv[(chain[0], chain[1], 2)]
                    elif self.verbose:
                        sys.stdout.write(
                            f"""Found model {self.ScFv[(chain[0], chain[1], 2)].model} chain {self.ScFv[(chain[0], chain[1], 2)].ident} V{self.ScFv[(chain[0], chain[1], 2)].type} region.\n"""
                        )
            self.Chains.update(self.ScFv)
            # Remove duplicates from the Chains dictionary
            for chain in self.ScFv:
                if (chain[0], chain[1]) in self.Chains:
                    del self.Chains[(chain[0], chain[1])]

        if self.verbose and self.ABchains:
            sys.stdout.write("Done...\n")

    def pair_chains(self):
        # Pair chains within the models found in the PDB file.
        # Use distance constraint between two highly conserved interface residues
        if self.verbose and self.ABchains:
            sys.stdout.write("Paring VH and VL...\n")
        Models = {}
        for chain in self.Chains:
            if self.Chains[chain].type in ["H", "L"]:
                try:
                    Models[chain[0]].append(chain)
                except KeyError:
                    Models[chain[0]] = [chain]

        # TODO implement function to compute distance between two chains
        self.Fvs = []
        for model in Models:
            for comb in itertools.combinations(Models[model], 2):
                if self.Chains[comb[0]].type == self.Chains[comb[1]].type:
                    continue
                # Use a distance cutoff of 20 Angstroms between two conserved interface atoms in order to pair domains.
                elif (
                    numpy.linalg.norm(
                        numpy.array(self.Chains[comb[0]].point)
                        - numpy.array(self.Chains[comb[1]].point)
                    )
                    < 20.0
                ):
                    self.Fvs.append(
                        Fv(self.Chains[comb[0]], self.Chains[comb[1]], self.name)
                    )
        if not self.Fvs:
            if self.ABchains:
                unpaired_chains = self.filepath, ", ".join([f'{c[1]} (V{self.Chains[c].type})' for c in self.ABchains])
                
                raise Exception(
                    f"""Could not locate Fvs in {self.filepath}
                    Unpaired antibody chains were found: {unpaired_chains}
                    These are either of all the same type (VH or VL) or do not meet the distance constraint.
                    Please use the scfv option for single chain fvs.\n"""
                )
            else:
                raise Exception(
                    f"""No antibody chains found in {self.filepath}.\n
                    check that heavy and light chains have chain identifiers H and L respectively.\n"""
                )

class Chain:
    """A class to describe a chain in a pdb file.
    Antibody numbering will be applied to each chain in the file. If it succeeds the coordinates will be annotated with chothia numbering."""

    def __init__(self, lines, cid, args):
        self.lines = lines
        self.model = cid[0]
        self.ident = cid[1]
        self.args = args
        self.annotate()

    def annotate(self):

        # Convert pdb lines to pir format. Returns the file handle.
        pirfile, Sequence, ResID = pdb2pir(self.lines)

        if len(Sequence) <= 2000:
            # Number using abnum either locally or online. OR let the user tell us that it is already chothia numbered.
            # The executable at AnnotationProgPath should take as input a .pir file and give the same format output as the abnum server to stdout.
            # Errors should go to stderr and the executable should be quiet if it cannot number the sequence given. i.e return '' to both stdout and stderr.
            # Only tested with Abysis
            if self.args.usernumbered:
                # Parse the PDB to get the correct annotation.
                ABnumOut = UserAnnotation(Sequence, ResID, self.ident)
            elif self.args.online or AnnotationProgPath == "Use online":
                ABnumOut = OnlineAnnotation(Sequence)
            elif AnnotationProgPath:
                subpr = subprocess.Popen(
                    [AnnotationProgPath, pirfile, "-c"],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                )
                ABnumOut = subpr.communicate()
            else:
                raise Exception(
                    "No numbering option specified and no numbering program found"
                )
        else:
            ABnumOut = [False, False]

        os.remove(pirfile)

        if ABnumOut[1]:
            raise Exception(ABnumOut[1])

        # Abnum is quiet if it fails - we assume this is a non-antibody chain - call it an antigen - but does not necesserily mean that it is!
        # Antigens are not defined in this way
        if not ABnumOut[0]:
            self.type = "Antigen"
            self.AGlines = self.lines
        else:
            ######################
            # Parse ABnum output #
            ######################
            try:
                Annotation = [
                    res
                    for res in map(str.split, ABnumOut[0].strip().split("\n"))
                    if res[1] != "-"
                ]
            except IndexError:
                raise Exception(
                    "Annotation failed. Unexpected annotation file format. Starts with: %s "
                    % ABnumOut[0][:50]
                )

            # Match the start of the Annotated sequence to the input sequence. Abnum leaves out insertions at the start of the sequence.
            # Map the annotation to the ResID. dictionary with residue id as key, annotation as values
            mapping, self.sequence = AlignSequence(Annotation, Sequence, ResID)

            # Assign the chain type: H, L
            self.type = Annotation[0][0][0]

            # Annotate the PDBfile with chothia numbering. Heavy and Light chain types will be labelled as H and L respectively
            self.FVATOMS, self.unannotated = [], []
            for line in self.lines:
                rid = line[22:27].strip()
                if rid in mapping:
                    try:
                        self.FVATOMS.append(
                            line[:21]
                            + self.type
                            + ("%d" % int(mapping[rid][1:])).rjust(4)
                            + " "
                            + line[27:]
                        )
                    except ValueError:
                        # If you have insertion ( e.g. H97A )
                        self.FVATOMS.append(
                            line[:21]
                            + self.type
                            + mapping[rid][1:-1].rjust(4)
                            + mapping[rid][-1]
                            + line[27:]
                        )
                # If we have started to find FVATOMS then add the rest to unannotated - important for scFvs
                else:  # if self.FVATOMS:
                    self.unannotated.append(line)
            # Extract an interface ca atom for pairing calculations.
            if self.type == "H":
                self.point = [
                    list(map(float, [l[30:38], l[38:46], l[46:54]]))
                    for l in self.FVATOMS
                    if l.split()[2] + self.type + l[22:27].strip() == "CAH37"
                ][0]
            elif self.type == "L":
                self.point = [
                    list(map(float, [l[30:38], l[38:46], l[46:54]]))
                    for l in self.FVATOMS
                    if l.split()[2] + self.type + l[22:27].strip() == "CAL87"
                ][0]


class Fv:
    """A class to compile chains which have been paired to form an Fv region"""

    def __init__(self, c1, c2, PDBname):
        self.c1 = c1
        self.c2 = c2
        self.sequence = {}
        self.sequence.update(c1.sequence)
        self.sequence.update(c2.sequence)
        self.GetName(PDBname)
        self.CompileLines()

    def GetName(self, PDBname):
        if self.c1.type == "H" and self.c2.type == "L":
            self.name = PDBname + "_" + self.c1.ident + self.c2.ident
            self.chains = {"H": self.c1, "L": self.c2}
        elif self.c1.type == "L" and self.c2.type == "H":
            self.name = PDBname + "_" + self.c2.ident + self.c1.ident
            self.chains = {"L": self.c1, "H": self.c2}
        else:
            raise Exception("Pairing of non-AB chains as an Fv")

        if self.c1.model:
            self.name += "_m%s" % self.c1.model

    def CompileLines(self):
        self.Lines = self.chains["H"].FVATOMS + self.chains["L"].FVATOMS


########################
# Annotation functions #
########################


def UserAnnotation(Sequence, ResID, ident):
    """The user gives a PDB containing an Fv region. It must be chothia numbered and be labelled with either H or L as the chain identifier"""
    annotation = ""
    if ident in ["H", "L"]:
        for i in range(len(Sequence)):
            annotation += ident + ResID[i] + " " + Sequence[i] + "\n"
    return annotation, ""


def OnlineAnnotation(sequence):
    """Use the public Abnum server for annotation.
    Note that this uses an PUBLIC webserver - this will only be done if the user requests it or explicitly sets it as their default method"""

    Url = (
        "http://www.bioinf.org.uk/cgi-bin/abnum/abnum.pl?plain=1&aaseq=%s&scheme=-c"
        % ("".join(sequence))
    )  # chothia numbering input
    Annfd, Annfile = tempfile.mkstemp(".dat", "annot")
    try:
        urllib.request.urlretrieve(Url, Annfile)
    except IOError:
        raise Exception(
            "Online annotation failed. Unable to access %s. Please check you have working internet connection and the ABnum server is functioning.\n"
            % Url
        )

    f = os.fdopen(Annfd, "r")
    annotation = f.read()
    f.close()
    if annotation == "\n":  # It failed
        annotation = ""
    os.remove(Annfile)
    return annotation, ""


def AlignSequence(annotated, target, resid):
    """Align the annotated sequence we get out of the annotation method and the sequence it was given. Needed so just in case there is a duff output from the annotation program."""
    n = len(annotated)
    N = len(target)
    tstore = "".join(target)
    MaxCut = N - n
    cut = 0
    while 1:
        aligned = True
        for i in range(n):
            if annotated[i][1] != target[i]:
                aligned = False
                target = target[1:]
                resid = resid[1:]
                cut += 1
                break
        if aligned:
            break
        if cut > MaxCut:
            raise Exception(
                "Could not align annotated sequence and sequence from PDB\n%s\n%s\n"
                % ("".join([annotated[i][1] for i in range(n)]), tstore)
            )
    mapping = dict((resid[i], annotated[i][0]) for i in range(n))
    sequence = dict((annotated[i][0], annotated[i][1]) for i in range(n))
    return mapping, sequence

def timeout(signum, frame):
    raise OSError

#########################
# Calculation functions #
#########################


def create_coreset(fname):
    """Parses the file so that it only contains the coreset residues in each domain"""
    try:
        # Parse the input file
        Hfd, Hf = tempfile.mkstemp(".pdb", "H")
        Lfd, Lf = tempfile.mkstemp(".pdb", "L")

        Htmp = os.fdopen(Hfd, "w")
        Ltmp = os.fdopen(Lfd, "w")
        fin = open(fname, "r").readlines()
        for line in fin:
            l = line.split()
            if not "ATOM" in line:
                continue
            elif l[4] == "L" and l[5] in coreset['light']:
                Ltmp.write(line)
            elif l[4] == "H" and l[5] in coreset['heavy']:
                Htmp.write(line)  
        Htmp.close()
        Ltmp.close()
    except Exception as exe:
        os.remove(Hf)
        os.remove(Lf)
        raise Exception(str(exe) + "\n")
    
    return Hf, Lf

def mapvectors(fname, PAPS_def=False):
    """Maps the reference frames (planes) onto to VH and VL domains of an Fv structure (fname is chothia numbered pdb file
    with VH as H chain and VL as L chain. PAPS_def means use the same definition of C that Abhinandan  and Martin did when calculating
    their torsion angle (makes HL should be the same as their packing angle as defined in authors' paper)"""
    # Get transformation matrices by aligning the core of the domains
    Hf, Lf = create_coreset(fname)
    
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

###########################################
# Parsing functions and control functions #
###########################################


def get_fvs(FileNames, Mode, args):
    """Function to extract all the Fv regions from a PDB file. It creates a new chothia numbered PDB file for each Fv region and returns the file handle for it"""
    FvRegions = {}
    StoreMode = {}
    FvSequences = {}

    for File in FileNames:
        
        Fvs = PDB(File, FileNames[File], args).Fvs
        for Fv in Fvs:
            if args.store:
                # Stores the FvRegion
                FvRegions[Fv.name] = os.path.join(
                    user_datapath, "user_fvs", Fv.name + ".pdb"
                )
                StoreMode[Fv.name] = Mode[File]
                FvSequences[Fv.name] = Fv.sequence
            else:
                # Stores the FvRegion in a tempory file.
                fid, FvRegions[Fv.name] = tempfile.mkstemp(".pdb", Fv.name)
                os.close(fid)
                FvSequences[Fv.name] = Fv.sequence
            Fvfile = open(FvRegions[Fv.name], "w")
            Fvfile.writelines(Fv.Lines)
            Fvfile.close()

    # Return the FvRegions. Dictionary with the name of the Fv as the key (name_HL) and the file handle to the Fv region pdb
    return FvRegions, StoreMode, FvSequences

def GetAngles(args):
    """Calculate the angles for all the Fv regions in the PDB structure files defined in the -i argument."""
    NewAngles = {}
    out = sys.stdout
    try:
        if args.store == "n":
            args.store = False
        elif not args.store:
            out.write("Store new data in userdata? (local storage only)\n")
            while 1:
                answer = input("Y/N: ")
                if answer.upper() in ["Y", "N"]:
                    break
            if answer.upper() == "Y":
                args.store = True

        if args.store:
            if not user_datapath:  # create a repository if it has not been created
                args.store = ChangeUserPath(args)  # allow for user to change mind.
            if args.store and not os.path.exists(
                os.path.join(user_datapath, "user_fvs")
            ):
                os.mkdir(os.path.join(user_datapath, "user_fvs"))

        if args.store:
            NewFileNames, StoreMode = CheckNamespace(args)
            NewFvRegions, StoreMode, FvSequences = get_fvs(NewFileNames, StoreMode, args)
        else:
            # FIXME 
            # NewFileNames = dict(
            #     zip(
            #         map(lambda x: os.path.splitext(os.path.split(x)[1])[0], args.i),
            #         args.i,
            #     )
            # )
            # temporary patch #
            NewFileNames = {args.name: args.i}
            ###################
            NewFvRegions, StoreMode, FvSequences = get_fvs(
                NewFileNames, {}, args
            )  # StoreMode will come back empty
        table = []
        if not args.q and NewFvRegions:
            table.append(["Structure", "HL", "HC1", "LC1", "HC2", "LC2", "dc"])

        PredictedAngles = {}
        for structure in NewFvRegions:

            # Calculate the angles and extract the sequence of the antibody.
            NewAngles[structure] = angles(NewFvRegions[structure])

            if args.store:
                AppendAngleStore(
                    structure,
                    NewAngles[structure],
                    FvSequences[structure],
                    StoreMode[structure],
                )
                AppendSequenceStore(
                    structure, FvSequences[structure], StoreMode[structure]
                )
            table.append(
                [structure]
                + [
                    NewAngles[structure][a]
                    for a in ["HL", "HC1", "LC1", "HC2", "LC2", "dc"]
                ]
            )

        # Print table to screen nicely
        if table:
            dataIO.print_table(table, out)

        if not args.store:
            for s in NewFvRegions:
                try:
                    os.remove(NewFvRegions[s])
                except OSError:
                    continue
        return NewAngles
    except Exception as exe:
        raise Exception("Angle calculation failed: " + str(exe) + "\n")

@dataclass
class GetAnglesargs:
    """class to simulate the arguments the GetAngles function would receive were it being called from the command line
    used in testing module output against benchmarks"""
    i: str 
    name: str
    scfv: List = None
    store: str = 'n'
    q: bool = True
    usernumbered: bool = True

def validate_angles(file, HC2, HC1, LC2, LC1, dc, HL, rel_tol=1e-2):
    name = file[:4]
    angle_keys = ['HC2', 'HC1', 'LC2', 'LC1', 'dc', 'HL']
    true_angles = [HC2, HC1, LC2, LC1, dc, HL]
    args = GetAnglesargs(examples/file, name)
    new_angles = GetAngles(args)[f'{name}_HL']
    assert all([math.isclose(new_angles[key], angle, rel_tol=rel_tol) for key, angle in zip(angle_keys, true_angles)])

if __name__ == '__main__':
    
    examples = data_path.parent/'data'/'example_pdbs'

    angles_4KQ3 = {'HC2':114.97, 'HC1':71.58, 'LC2':83.15, 'LC1':119.49, 'dc':16.00, 'HL':-61.10}
    validate_angles('4KQ3_Fv.pdb', **angles_4KQ3)

    angles_2ATK = {'HL': -54.09, 'HC1': 70.76, 'HC2': 114.09, 'LC1': 123.29, 'LC2': 83.42, 'dc': 16.48}
    validate_angles('2ATK_Fv.pdb', **angles_2ATK)

    angles_1U8L = {'HL': -56.41, 'HC1': 68.73, 'HC2': 118.87, 'LC1': 123.11, 'LC2': 82.99, 'dc': 15.88}
    validate_angles('1U8L_Fv.pdb', **angles_1U8L)
