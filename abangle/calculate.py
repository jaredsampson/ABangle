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

# Get coreset positions
coresetL = [
    l.strip()[1:] for l in open(os.path.join(data_path, "Lcoresetfw.txt")).readlines()
]
coresetH = [
    l.strip()[1:] for l in open(os.path.join(data_path, "Hcoresetfw.txt")).readlines()
]

# Read in the plane vectors precalculated on the consensus structure
Lpos = [
    [float(num) for num in line.split()] 
    for line in (data_path/'pcL.txt').open()
]

Hpos = [
    [float(num) for num in line.split()] 
    for line in (data_path/'pcH.txt').open()
]

# Amino acid translation
AA = [
    "R",
    "H",
    "K",
    "D",
    "E",
    "S",
    "T",
    "N",
    "Q",
    "C",
    "G",
    "P",
    "A",
    "V",
    "I",
    "L",
    "M",
    "F",
    "Y",
    "W",
]
aa3 = [
    "ARG",
    "HIS",
    "LYS",
    "ASP",
    "GLU",
    "SER",
    "THR",
    "ASN",
    "GLN",
    "CYS",
    "GLY",
    "PRO",
    "ALA",
    "VAL",
    "ILE",
    "LEU",
    "MET",
    "PHE",
    "TYR",
    "TRP",
]
aa3to1 = dict(list(zip(aa3, AA)))

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


def pdb2pir(lines):
    """Create a .pir file which is used by a local abnum - available under licence of Abysis. Extract the sequence and the residue identities"""
    pirfd, pirfile = tempfile.mkstemp(".pir")
    f = os.fdopen(pirfd, "w")
    seq = []
    ResID = []
    for line in lines:
        l = line.split()
        if line[13:16].strip() != "CA":
            continue
        resid = line[22:27].strip()
        if resid in ResID:
            # Multiple occupancy - ignore any but the first occurance
            continue
        aa = line[16:20].strip()
        try:
            seq.append(aa3to1[aa])  # try to append the one letter code to the sequence.
        except KeyError:  # if we fail then there may be multiple occupancy
            if len(aa) == 4:  # check for residue name such as ASER or BSER
                aa3 = aa[1:]
                try:
                    seq.append(aa3to1[aa3])
                except KeyError:
                    sys.stderr.write(
                        "Warning: skipping atom: Unknown aa type: %s\n" % aa
                    )
                    continue
            else:
                sys.stderr.write("Warning: skipping atom: Unknown aa type: %s\n" % aa)
                continue
        ResID.append(resid)
    # Output to pir file.
    f.write(">P1;abchain\n")
    f.write("pir file generated by ABangle.calculate.pdb2pir\n")
    for aa in seq:
        f.write(aa)
    f.write("*\n")
    f.close()
    # Return the pir file handle, the sequence extracted and the residue IDs which correspond to the sequence
    return pirfile, seq, ResID


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
            elif l[4] == "L" and l[5] in coresetL:
                Ltmp.write(line)
            elif l[4] == "H" and l[5] in coresetH:
                Htmp.write(line)  
        Htmp.close()
        Ltmp.close()
    except Exception as exe:
        os.remove(Hf)
        os.remove(Lf)
        raise Exception(str(exe) + "\n")

    with open(Hf) as f:
        with open('Hf.pdb', 'w') as w:
            w.write(f.read())
    
    return Hf, Lf


def mapvectors(fname, PAPS_def=False):
    """Maps the reference frames (planes) onto to VH and VL domains of an Fv structure (fname is chothia numbered pdb file
    with VH as H chain and VL as L chain. PAPS_def means use the same definition of C that Abhinandan  and Martin did when calculating
    their torsion angle (makes HL should be the same as their packing angle as defined in authors' paper)"""
    # Get transformation matrices by aligning the core of the domains
    Hf, Lf = create_coreset(fname)
    
    uL = align(os.path.join(data_path, "consensus_L.pdb"), Lf)
    uH = align(os.path.join(data_path, "consensus_H.pdb"), Hf)
    # os.remove(Hf)
    # os.remove(Lf)

    if PAPS_def:
        # The centroids of interface residues.
        cH = Hpos[2]
        cL = Lpos[2]
    else:
        # The minimally varying centroid vector is at. As calculated.:
        cH = [
            -10 * 0.5 * Hpos[0][i] + 1 * 0.5 * Hpos[1][i] + Hpos[2][i] for i in range(3)
        ]
        cL = [
            6 * 0.5 * Lpos[0][i] - 2 * 0.5 * Lpos[1][i] + Lpos[2][i] for i in range(3)
        ]

    # Define the plane vectors from the centroid point
    # On VL domain
    L1 = [cL[i] + Lpos[0][i] for i in range(3)]
    L2 = [cL[i] + Lpos[1][i] for i in range(3)]

    # On VH domain
    H1 = [cH[i] + Hpos[0][i] for i in range(3)]
    H2 = [cH[i] + Hpos[1][i] for i in range(3)]

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


def H3LoopLength(seq):
    # Chothia definition of H3 - residues and insertions
    H3 = ["H95", "H96", "H97", "H98", "H99", "H100", "H101", "H102"]
    n = 0
    for res in seq:
        try:
            int(res[-1])
            if res in H3:
                if seq[res] != "-":
                    n += 1
        except ValueError:
            # The residue numbering ends in a letter.
            if res[:-1] in H3:
                if seq[res] != "-":
                    n += 1
    return n


def L1LoopLength(seq):
    # Chothia definition of L1 - residues and insertions
    L1 = ["L24", "L25", "L26", "L27", "L28", "L29", "L30", "L31", "L32", "L33", "L34"]
    n = 0
    for res in seq:
        try:
            int(res[-1])
            if res in L1:
                if seq[res] != "-":
                    n += 1
        except ValueError:
            # The residue numbering ends in a letter.
            if res[:-1] in L1:
                if seq[res] != "-":
                    n += 1
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


def ChangeUserPath(args):
    """Function to change or create the user repository path.  This is where all the user's data is
    stored."""
    global user_datapath

    if user_datapath:
        sys.stdout.write("Current user_datapath is: %s\n" % user_datapath)
    elif savedpath:
        sys.stdout.write("Saved user_datapath: %s was not found\n" % savedpath)
    p = input("Please provide a path to place a user repository (s to skip):\n")
    if p.lower() == "s":
        return False

    newpath = os.path.abspath(os.path.join(p, "ABangleData/"))
    # Handle some potential errors - this may not be completely robust.
    if not os.path.exists(newpath):
        try:
            os.mkdir(newpath)
        except OSError as exe:
            if str(exe).startswith("[Errno 13]"):
                sys.stderr.write("No write privelages for %s.\n" % os.path.abspath(p))
            else:
                sys.stderr.write(
                    "Path %s does not exist. Please provide an existing path to create a repository\n"
                    % os.path.abspath(p)
                )
            return False
    elif not (os.access(newpath, os.R_OK) and os.access(newpath, os.W_OK)):
        sys.stderr.write("No read/write privelages for %s.\n" % newpath)
        return False

    if not os.path.exists(os.path.join(newpath, "user_fvs")):
        try:
            os.mkdir(os.path.join(newpath, "user_fvs"))
        except OSError as exe:
            if str(exe).startswith("[Errno 13]"):
                sys.stderr.write("No write privelages for %s.\n" % os.path.abspath(p))
            return False
    elif not (os.access(newpath, os.R_OK) and os.access(newpath, os.W_OK)):
        sys.stderr.write("No read/write privelages for %s.\n" % newpath)
        return False

    user_datapath = newpath
    ufname = open(os.path.join(path, "config/userdatapath.txt"), "w")
    ufname.write(user_datapath)
    ufname.close()
    # Create the data store files.
    CreateStore()

    return True


def ChangeNumberingProgPath(args):
    global AnnotationProgPath

    if AnnotationProgPath:
        sys.stdout.write(
            "Current path to antibody numbering script is: %s\n" % AnnotationProgPath
        )

    while 1:
        a = input(
            "Use public webserver numbering (http://www.bioinf.org.uk/abs/abnum/) for all numbering? y/n/q\n"
        )
        if a.lower() == "q":
            raise Exception("Quit\n")
        elif a.lower() in ["y", "n"]:
            break
    if a.lower() == "y":
        newpath = "Use online"
    else:
        while 1:
            sys.stdout.write(
                "Only abysis is currently supported for local numbering.\n"
            )
            p = input(
                "Provide path to numbering script (kabnum_wrapper.pl for abysis)? y/n/q\n"
            )
            if p.lower() == "q":
                raise Exception("Quit\n")
            elif p.lower() in ["y", "n"]:
                break
        if p.lower() == "n":
            raise Exception(
                "User must either allow online numbering or local numbering with abysis. Quitting\n"
            )
        else:
            newpath = input(
                "Path to numbering script (kabnum_wrapper.pl for abysis): "
            )
            if newpath.endswith("kabnum_wrapper.pl"):
                newpath = os.path.split(newpath)[0]

            newpath = os.path.abspath(os.path.join(newpath, "kabnum_wrapper.pl"))
            if not os.path.exists(newpath):
                raise Exception("%s was not found\n" % newpath)

    AnnotationProgPath = newpath
    Progfname = open(os.path.join(path, "config/AnnotationProgPath.txt"), "w")
    Progfname.write(AnnotationProgPath)
    Progfname.close()


def DeleteStore(remove):
    """Function to delete entries for the datafiles"""
    if not user_datapath:
        raise Exception("No user data found")

    if remove[0] == "all":
        remove = list(dataIO.load(
            os.path.join(user_datapath, "UserAngles.dat"), header=True, rownames=0
        )[0].keys())
        # This wipes the userdata angle store and rewrites the header
        CreateStore()
        sys.stdout.write("User's data deleted\n")
        # Delete the files in the USERDATA/USERFVS directory.
    else:
        # remove the structures requested one by one. If the exact file is not requested then it will not be removed. Eg. -clear AB1 will not remove AB1_HL and AB1_BC but -clear AB1_HL AB1_BC will do.
        for s in remove:
            AppendAngleStore(s, {}, "", "r", remove=True)
            AppendSequenceStore(s, {}, "r", remove=True)

    for s in remove:
        fname = "%s.pdb" % s
        pathname = os.path.abspath(os.path.join(user_datapath, "user_fvs", fname))
        try:
            os.remove(pathname)
        except:
            sys.stderr.write(
                "Warning: Structure %s was not found in the users' store\n" % fname
            )


def CreateStore():
    """Create storage files for user data. This stores the orientation measures and the angles for the data that the user has inputted (and saved/stored).
    These files are in the same format as Angles.dat and Sequences.dat in the abangle/data/ directory."""
    Out = open(os.path.join(user_datapath, "UserAngles.dat"), "w")
    # 	Out.write( "Code\tHL\tHC1\tLC1\tHC2\tLC2\tdc\tMethod\tRes\tRfac\tBstate\tFvCont\tAGtype\tAGsize\tHgroup\tLgroup\tSpecies\tLtype\tH3length\tL1length\n" )
    Out.write(
        "Code\tHL\tHC1\tLC1\tHC2\tLC2\tdc\tMethod\tRes\tRfac\tBstate\tAGtype\tHgroup\tLgroup\tSpecies\tLtype\tH3length\tL1length\n"
    )
    Out.close()
    Out = open(os.path.join(user_datapath, "UserSequences.dat"), "w")
    Out.writelines(["Code"] + ["\t%s" % res for res in Residues] + ["\n"])
    Out.close()


def AppendAngleStore(structure, dat, seq, mode, remove=False):
    """Add (or remove) a user structures's angles."""
    try:
        UserAngles = open(os.path.join(user_datapath, "UserAngles.dat"), mode)
    except IOError:
        # File does not exist. First run or user has deleted manually. Create a new one and open again.
        CreateStore()
        UserAngles = open(os.path.join(user_datapath, "UserAngles.dat"), mode)

    if not remove:
        dat["H3length"] = H3LoopLength(seq)
        dat["L1length"] = L1LoopLength(seq)
    try:
        lines = UserAngles.readlines()
        # opened in read mode - overwrite record for stored data.
        # 		outlines=["Code\tHL\tHC1\tLC1\tHC2\tLC2\tdc\tMethod\tRes\tRfac\tBstate\tFvCont\tAGtype\tAGsize\tHgroup\tLgroup\tSpecies\tLtype\tH3length\tL1length\n"]
        outlines = [
            "Code\tHL\tHC1\tLC1\tHC2\tLC2\tdc\tMethod\tRes\tRfac\tBstate\tAGtype\tHgroup\tLgroup\tSpecies\tLtype\tH3length\tL1length\n"
        ]

        for line in lines:
            l = line.split()[0]
            if l == "Code":
                continue
            elif l == structure:
                continue
            else:
                outlines.append(line)
        if not remove:
            outlines.append(
                structure
                + "\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t?\t?\t?\t?\t?\t?\t?\t?\t?\t%d\t%d\n"
                % tuple(
                    dat[a]
                    for a in [
                        "HL",
                        "HC1",
                        "LC1",
                        "HC2",
                        "LC2",
                        "dc",
                        "H3length",
                        "L1length",
                    ]
                )
            )
        UserAngles.close()

        # Re-open in write mode and write out all lines
        Out = open(os.path.join(user_datapath, "UserAngles.dat"), "w")
        Out.writelines(outlines)
        Out.close()

    except IOError:
        # Opened in append mode
        UserAngles.write(
            structure
            + "\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t?\t?\t?\t?\t?\t?\t?\t?\t?\t%d\t%d\n"
            % tuple(
                dat[a]
                for a in [
                    "HL",
                    "HC1",
                    "LC1",
                    "HC2",
                    "LC2",
                    "dc",
                    "H3length",
                    "L1length",
                ]
            )
        )
        UserAngles.close()


def AppendSequenceStore(structure, seq, mode, remove=False):
    """Add (or remove) a user structures's sequence. The storage method may be slow if doing many structures at a time."""
    try:
        # Grap the UserSequences file.
        UserSequencesDict, UserResidues = dataIO.load(
            os.path.join(user_datapath, "UserSequences.dat"),
            header=True,
            rownames=0,
            conv=False,
        )
    except IOError:
        CreateStore()
        UserSequencesDict, UserResidues = {}, Residues

    # Check that there is a column for each of the residues that have been annotated in the new sequence. Change the mode to 'r' if additional column is required.
    for r in list(seq.keys()):
        if r not in UserResidues:
            UserResidues.append(r)
            mode = "r"

    # Open the file in append mode or read mode.
    UserSequencesFile = open(os.path.join(user_datapath, "UserSequences.dat"), mode)

    try:
        lines = UserSequencesFile.readlines()
        # opened in read mode - overwrite record for stored data.
        outlines = ["Code\t" + "\t".join(res for res in UserResidues) + "\n"]
        for line in lines:
            code = line.split()[0]
            if code == "Code":
                continue
            elif code == structure:
                continue
            else:
                outlines += [
                    code
                    + "\t"
                    + "\t".join(
                        UserSequencesDict[code][res]
                        if res in UserSequencesDict[code]
                        else "-"
                        for res in UserResidues
                    )
                    + "\n"
                ]
        if not remove:
            outlines += [
                structure
                + "\t"
                + "\t".join(seq[res] if res in seq else "-" for res in UserResidues)
                + "\n"
            ]
        UserSequencesFile.close()

        # Re-open in write mode and write out all lines
        Out = open(os.path.join(user_datapath, "UserSequences.dat"), "w")
        Out.writelines(outlines)
        Out.close()

    except IOError:
        # Opened in append mode
        UserSequencesFile.write(
            structure
            + "\t"
            + "\t".join(seq[res] if res in seq else "-" for res in UserResidues)
            + "\n"
        )
        UserSequencesFile.close()


def CheckNamespace(args):
    """Check for namespace clashes in the stored data"""
    # Load the current NameSpace of the user's structures
    try:
        CurrentNames = []
        for name in [
            n.split("_")
            for n in list(dataIO.load(
                os.path.join(user_datapath, "UserAngles.dat"), header=True, rownames=0
            )[0].keys())
        ]:
            if len(name) == 2:
                CurrentNames.append(name[0])
            elif name[-1][0] == "m":
                CurrentNames.append("_".join(name[:-2]))
    except IOError:
        raise Exception(
            "%s was not found. Use -change_data_path to create a new repository"
            % (os.path.join(user_datapath, "UserAngles.dat"))
        )

    # Structures with multiple models will be stored as name_HL_mi. Where name is the PDB file name, HL the heavy and light chain identifiers and i the model number.
    # Structures with no model separators in the file will be stored as name_HL.

    # Check whether we would have a NameSpace clash. Offer to overwrite, rename  or skip the structure.
    NewFileNames = {}
    StoreMode = {}

    PathsFiles = list(map(os.path.split, args.i))
    names = [os.path.splitext(x[1])[0] for x in PathsFiles]

    for n in range(len(names)):
        if names[n] in CurrentNames:
            sys.stdout.write(
                "The name %s already exists. Overwrite (o), rename (r) or skip structure (s) ?\n"
                % names[n]
            )
            while 1:
                answer = input("o/r/s: ").lower()
                if answer in "ors":
                    break
            sys.stdout.write("\n")
            if answer == "o":
                NewFileNames[names[n]] = args.i[n]
                StoreMode[names[n]] = "r"
            elif answer == "r":
                while 1:
                    newname = input(
                        "\nNew name (no file extension). Leave blank to skip structure: "
                    )
                    if not newname:
                        answer = "s"
                        break
                    elif newname in CurrentNames or newname in names:
                        sys.stdout.write(
                            "\nName %s already exists or in -i/f argument. " % newname
                        )
                    else:
                        NewFileNames[newname] = args.i[n]
                        StoreMode[newname] = "a"
                        break
            if answer == "s":
                break
        else:
            # If does not exist, then let through as is.
            NewFileNames[names[n]] = args.i[n]
            StoreMode[names[n]] = "a"

    # return the dictionary with New Namespace with values as file names
    return NewFileNames, StoreMode


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
            # NewFileNames = dict(
            #     zip(
            #         map(lambda x: os.path.splitext(os.path.split(x)[1])[0], args.i),
            #         args.i,
            #     )
            # )
            NewFileNames = {args.name: args.i}
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
    
    examples = data_path/'example_pdbs'

    angles_4KQ3 = {'HC2':114.97, 'HC1':71.58, 'LC2':83.15, 'LC1':119.49, 'dc':16.00, 'HL':-61.10}

    validate_angles('4KQ3_abnum.pdb', **angles_4KQ3)
    

