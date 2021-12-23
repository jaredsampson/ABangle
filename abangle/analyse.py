"""
DESCRIPTION	
	
	abangle.analyse.py
	A module to analyse the orientation between the VH and VL domains in an Antibody Fv regions.

OUTPUT
	
	~ Individual or sets of structures and their orientation measures.
	~ Text information to stdout
	~ Plots of the orientation measures for requested structures
	~ Interface with PyMOL to visualise structures.

AUTHOR
	
	2012
	James Dunbar - Oxford protein informatics group.
	james.dunbar@dtc.ox.ac.uk
	www.stats.ox.ac.uk/~dunbar/
		
	Supervisors
	Prof C.Deane - Oxford protein informatics group.
	Dr Angelika Fuchs (Roche) and Dr Jiye Shi (UCB Celltech)\n
"""

#########################
# Import python modules #
#########################

from operator import mul
from math import log
from textwrap import wrap
import argparse, sys, random, os, subprocess, tempfile, re
import itertools

##########################
# Import abangle modules #
##########################

from abangle import align, dataIO

####################
# Global variables #
####################

# Load the data files
path = os.path.split(__file__)[0]
data_path = os.path.join(path, "data")

# This is the redundant set. Structures can be requested from it
All_Angles, anglenames = dataIO.load(
    os.path.join(data_path, "All_Angles.dat"), header=True, rownames=0
)
All_Sequences, Residues = dataIO.load(
    os.path.join(data_path, "All_sequences.dat"), header=True, rownames=0, conv=False
)
All_Names = (
    {}
)  # parse all the names with the same root name ( e.g. 12E8_HL and 12E8_PM have 12E8 as the root name )
for name in All_Angles.keys():
    try:
        All_Names[name.split("_")[0]].append(name)
    except KeyError:
        All_Names[name.split("_")[0]] = [name]


# This is the non-redundant set. This is used to search for structures and plot the background by default.
Sequences, Residues = dataIO.load(
    os.path.join(data_path, "Sequences.dat"), header=True, rownames=0, conv=False
)
Angles, anglenames = dataIO.load(
    os.path.join(data_path, "Angles.dat"), header=True, rownames=0
)

# This is the user's angles and sequences. It is updated as the user adds more structures and stores them.
try:
    user_datapath = (
        open(os.path.join(path, "config/userdatapath.txt")).readline().strip()
    )
    User_Sequences, userResidues = dataIO.load(
        os.path.join(user_datapath, "UserSequences.dat"),
        header=True,
        rownames=0,
        conv=False,
    )  # This will update as necessary.
    User_Angles, useranglenames = dataIO.load(
        os.path.join(user_datapath, "UserAngles.dat"), header=True, rownames=0
    )
    Users_Names = (
        {}
    )  # parse the names of the structures the user has saved so that they can use the root name to recall them.
    for name in [n.split("_") for n in User_Angles.keys()]:
        if len(name) == 2:
            try:
                Users_Names[name[0]].append("_".join(name))
            except KeyError:
                Users_Names[name[0]] = ["_".join(name)]
        elif name[-1][0] == "m":
            try:
                Users_Names["_".join(name[:-2])].append("_".join(name))
            except KeyError:
                Users_Names["_".join(name[:-2])] = ["_".join(name)]
except IOError:
    user_datapath = ""
    User_Sequences, userResidues = {}, []
    User_Angles, useranglenames = {}, []
    Users_Names = {}


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
    "-",
]
aa3 = [
    "arg",
    "his",
    "lys",
    "asp",
    "glu",
    "ser",
    "thr",
    "asn",
    "gln",
    "cys",
    "gly",
    "pro",
    "ala",
    "val",
    "ile",
    "leu",
    "met",
    "phe",
    "tyr",
    "trp",
    "-",
]
aadict = dict([(AA[i], aa3[i]) for i in range(len(AA))])  # for pymol.


class StructureSet:
    """A class for defining a set of structures requested on the command line. Stores the namespace and searches for the strucutrues requested."""

    def __init__(self, r, wa, f, wp, structures=[], label=""):
        self.factortype = {
            "Method": "c",
            "Bstate": "c",
            "AGtype": "c",
            "Hgroup": "c",
            "Lgroup": "c",
            "Ltype": "c",
            "Species": "c",
            "Res": "n",
            "Rfac": "n",
            "L1length": "n",
            "H3length": "n",
            "HL": "n",
            "HC1": "n",
            "HC2": "n",
            "LC1": "n",
            "LC2": "n",
            "dc": "n",
        }
        # Define the query parameters for this set.
        self.r = r
        self.wa = wa  # list of "with amino acids" - or's are not parsed yet
        self.f = f  # list of factors
        self.Classf = []  # Classifiable factors (e.g. species)
        self.Numf = []  # Numeric factors (e.g. angles)
        self.wp = wp  # list of "with factors" - or's are not parsed yet
        self.color = ""
        self.plotter = None
        # Hierarchy of things to search for might need to change for average speed.
        if not structures:
            self.SearchFactors()
            self.SearchResidues()
            self.GetLabel()
        else:
            self.structures = structures
            self.GetAngles()

            if label:
                self.Label = "User defined: " + "\n".join(wrap(label, width=50))
            elif len(self.structures) > 20:
                self.Label = "User defined: " + "\n".join(
                    wrap(" ".join(self.structures[:20] + [" ...."]), width=50)
                )
            else:
                self.Label = "User defined: " + "\n".join(
                    wrap(" ".join(self.structures), width=50)
                )

    def __str__(self):
        string = "Structure(s) with:\n"
        string += " ".join(self.Label.split("\n")) + "\n"
        if self.structures:
            for struc in self.structures:
                string += "%s\n" % struc
        else:
            string += "None"
        string += "\n"
        return string

    def __len__(self):
        return len(self.structures)

    def SearchFactors(self):
        """Search the structures to find those with a certain combination of factors. Split these into classifiable factors (e.g. Species, Hgroup etc) and numeric factors
        (e.g all the angles and dc as well as rfactor etc"""
        if self.wp:
            # Search for structures which fit the class parameters
            self.Classf = [f for f in self.f if self.factortype[f] == "c"]
            self.Classwp = Parselogic(
                [
                    self.wp[i]
                    for i in range(len(self.wp))
                    if self.factortype[self.f[i]] == "c"
                ]
            )

            if self.Classwp[0]:
                names = [
                    name
                    for name in Angles.keys()
                    if tuple(
                        Angles[name][f] for f in self.f if self.factortype[f] == "c"
                    )
                    in self.Classwp
                ]
            else:
                names = Angles.keys()

            # Filter these with those which meet the a numerical parameters.
            self.Numf = [f for f in self.f if self.factortype[f] == "n"]
            self.Numwp = Parselogic(
                [
                    self.wp[i]
                    for i in range(len(self.wp))
                    if self.factortype[self.f[i]] == "n"
                ]
            )
            if self.Numwp[0]:
                self.structures = []
                for searchfor in self.Numwp:
                    X, dx, funcs = InterpretNumeric(
                        [
                            ":".join([self.Numf[i], searchfor[i]])
                            for i in range(len(self.Numf))
                        ]
                    )

                    for name in names:
                        # Only take structures which satisfy all the conditions.
                        if reduce(
                            mul,
                            (funcs[f](Angles[name][f], X[f], dx[f]) for f in X.keys()),
                        ):
                            self.structures.append(name)
            else:
                self.structures = names
        else:
            self.structures = Angles.keys()

    def SearchResidues(self):
        """Search the structures to find those with a certain combination of amino acids at given residue positions"""
        if self.wa:

            self.Residueswa = Parselogic(self.wa)
            N = []
            for searchfor in self.Residueswa:
                for s in self.structures:
                    if tuple(Sequences[s][r] for r in self.r) == searchfor:
                        N.append(s)
            self.structures = N

    def GetAngles(self):
        """Retrieve the angles of the structures requested. If the user has been lazy with defining the exact Fv, retrieve all the Fvs within the original pdb requested.
        e.g 12E8 would pull out 12E8_HL and 12E8_PM"""
        self.Angles = dict((Ang, []) for Ang in anglenames)
        structures = (
            []
        )  # If the user does not include the chain identifier _HL then search for the root name.

        # Search the users angles first. Then the default.
        # Restructuring this.
        # Look up exact name in user angles
        # Look up exact name in pre-calculated angles
        # Look for root name in user angles
        # Look up root name in pre-calculated angles

        for code in self.structures:
            matches = []
            # look up the full name in user angles
            try:
                matches.append(User_Angles[code])
                structures.append(code)
            except KeyError:
                # look up the full name in pre-calculated structures
                try:
                    matches.append(All_Angles[code.upper()])
                    structures.append(code.upper())
                except KeyError:
                    # look up the the root name in User_Names and return all matching names - one loop at top.
                    try:
                        for name in Users_Names[code]:
                            matches.append(User_Angles[name])
                            structures.append(name)
                    except KeyError:
                        # look up the root name in All_Names and return all matching names
                        try:
                            for name in All_Names[code.upper()]:
                                matches.append(All_Angles[name])
                                structures.append(name)
                        except KeyError:
                            raise Exception("Structure %s not known\n" % code)

            for Ang in anglenames:
                for match in matches:
                    self.Angles[Ang].append(match[Ang])
        self.structures = structures

    def GetLabel(self):
        """Interpret the options for the structure set and produce a string to label it on output"""
        self.Label = ""
        if self.r:
            seen = dict((r, []) for r in self.r)
            for r in range(len(self.r)):
                substrings = []
                for a in self.Residueswa:
                    if a[r] not in seen[self.r[r]]:
                        substrings.append("%s " % a[r])
                        seen[self.r[r]].append(a[r])
                self.Label += "or ".join(substrings)
                self.Label += "at %s; " % self.r[r]
            if self.Classf or self.Numf:
                self.Label += "and "
        if self.Classf:
            seen = dict((f, []) for f in self.Classf)
            for f in range(len(self.Classf)):
                self.Label += "%s of " % self.Classf[f]
                substrings = []
                for p in self.Classwp:
                    if p[f] not in seen[self.Classf[f]]:
                        substrings.append("%s " % p[f])
                        seen[self.Classf[f]].append(p[f])
                self.Label += "or ".join(substrings)
            if self.Numf:
                self.Label += "and "
        if self.Numf:
            seen = dict((f, []) for f in self.Numf)
            for f in range(len(self.Numf)):
                self.Label += "%s " % self.Numf[f]
                substrings = []
                for p in self.Numwp:
                    if p[f] not in seen[self.Numf[f]]:
                        np = p[f].split(":")
                        if np[0] == "gt":
                            substrings.append("> %s " % str(np[1]))
                        elif np[0] == "lt":
                            substrings.append("< %s " % str(np[1]))
                        else:
                            substrings.append(
                                "+/- %s of %s " % (str(np[0]), str(np[1]))
                            )
                        seen[self.Numf[f]].append(p[f])
                self.Label += "or ".join(substrings)

        if not self.Label:
            self.Label = "all"
        self.Label += " : %d structures" % len(self.structures)
        # Split the label over 50 character lines for plotting legend
        self.Label = "\n".join(
            ["\n".join(wrap(block, width=50)) for block in self.Label.splitlines()]
        )


###################################
# Option interpretation functions #
###################################


def useall():
    """Modifies the global variables so that all the data is used. (Not just non-redundant set)"""
    global Angles
    global Sequences
    global Residues
    Angles = All_Angles
    Sequences, Residues = dataIO.load(
        os.path.join(data_path, "All_sequences.dat"),
        header=True,
        rownames=0,
        conv=False,
    )


def usebg(args):
    global Angles
    global Sequences
    Angles, anglenames = dataIO.load(args.bg, header=True, rownames=0)
    try:
        Sequences = dict((struc, All_Sequences[struc]) for struc in Angles)
    except:
        pass


def iswithin(angle, X, dx):
    """iswithin function. Is angle within dx of X?"""
    try:
        if abs(float(angle) - X) <= dx:
            return True
        else:
            return False
    except ValueError:  # Handle if we have missing data in numerical column
        return False


def isgreater(angle, X, dx):
    """isgreater function. Is angle greater than X?"""
    try:
        if float(angle) > X:
            return True
        else:
            return False
    except ValueError:
        return False


def isless(angle, X, dx):
    """isless function. Is angle less than X?"""
    try:
        if float(angle) < X:
            return True
        else:
            return False
    except ValueError:
        return False


def isequal(angle, X, dx):
    """isequal function. Is angle equal X?"""
    try:
        if float(angle) == X:
            return True
        else:
            return False
    except ValueError:
        return False


def InterpretNumeric(L):
    """Interpret a list of numeric factors"""
    arguments = [i.split(":") for i in L]
    X = dict((arg[0], float(arg[2])) for arg in arguments)
    func = {}
    dx = {}
    for arg in arguments:
        try:
            dx[arg[0]] = float(arg[1])
            func[arg[0]] = iswithin
        except ValueError:
            if arg[1] == "gt":
                func[arg[0]] = isgreater
            elif arg[1] == "lt":
                func[arg[0]] = isless
            elif arg[1] == "eq":
                func[arg[0]] = isequal
            else:
                raise Exception(
                    'Error: Unrecognised comparator or threshold "%s"' % arg[1]
                )
            dx[arg[0]] = None
    return X, dx, func


def ParseW(rawarg):
    """Parse the -wa or -wp arguments (with arguments) to give nested lists of with arguments"""
    return map(lambda x: x.split(), rawarg.split(","))


def Parselogic(l):
    """Parse a sequence of query parameters: + => 'or'.
    If you put ~ infront then parse to give 'not'"""
    w = map(lambda x: x.split("+"), l)
    wparsed = []
    for r in w:
        notFlag = [x.startswith("not") for x in r]
        if not sum(notFlag):
            wparsed.append(r)
            continue
        else:
            notAA = [r[i].strip("not") for i in range(len(r)) if notFlag[i]]

            wparsed.append([aa for aa in AA if aa not in notAA])
    return list(itertools.product(*wparsed))  #  find all combinations


####################
# Dependency check #
####################


def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


#######################
# Searching functions #
#######################
def target(args):
    """Find the structures which lie within the thresholds for structural similarity found using unbound homologues"""
    # Check that the user has not put in any angular aguments in the query.
    if set(args.f) & set(["HL", "HC1", "LC1", "HC2", "LC2"]):
        raise Exception("Angular parameters may not be used with -target option\n")

    if User_Angles.has_key(args.target):
        TargetAngles = User_Angles[args.target]
    elif All_Angles.has_key(args.target):
        TargetAngles = All_Angles[args.target]
    else:
        raise Exception("Structure %s not known\n" % args.target)

    args.f += ["HL", "HC1", "LC1", "HC2", "LC2"]
    # 	args.f += ['HC1', 'LC1', 'HC2', 'LC2']
    if args.wp:
        wp = args.wp.split(",")
    else:
        wp = [""]

    for i in range(len(wp)):
        wp[i] += " 4:%.2f 1.25:%.2f 1.25:%.2f 1.75:%.2f 2.25:%.2f" % tuple(
            TargetAngles[ang] for ang in ["HL", "HC1", "LC1", "HC2", "LC2"]
        )
    # 		wp[i] += " 1.25:%.2f 1.25:%.2f 1.75:%.2f 2.25:%.2f"%tuple( TargetAngles[ang] for ang in ['HC1', 'LC1', 'HC2', 'LC2'] )
    args.wp = ",".join(wp)

    return args


def GetFactorCombinations(Setr, Setwa, Setf, Setwp, args):
    """Find the factor combinations. These may be either amino acids at residue postions or factors.
    Returns a dictionary with the parameter combination and its observed rate and a dictionary with the combination's expected rate."""
    nr, nwa, nf, nwp = len(Setr), len(Setwa), len(Setf), len(Setwp)

    # Find out if there is missing arguments. If so find the frequency of missing things.
    ResNames, FacNames = [], []

    if Setr and (nr > nwa):
        ResNames = Setr[nwa:]
        Setr = Setr[:nwa]
        Setwa = Setwa

    if Setf and (nf > nwp):
        FacNames = Setf[nwp:]
        Setf = Setf[:nwp]
        Setwp = Setwp

    # Parse which structures to look analyse
    # If the user has supplied an arbritary set of structures in args.s options then use those.
    if args.s:
        Set = UserSets(args.s)[0]
    else:
        Set = StructureSet(Setr, Setwa, Setf, Setwp)

    Names = Set.structures

    # Find the property frequencies at each of the factor types and the combinations present.
    rf = dict((r, dict((a, 0) for a in AA)) for r in ResNames)
    ff = dict((f, {}) for f in FacNames)
    combinations = {}
    for struc in Names:
        try:
            comb = tuple(
                [All_Sequences[struc][r] for r in ResNames]
                + [All_Angles[struc][f] for f in FacNames]
            )
        except KeyError, e:
            if str(e).strip("'") != struc:
                raise Exception("Unrecognised factor %s" % e)

            # Struc might have come from the user's own set so look in them
            try:
                comb = tuple(
                    [User_Sequences[struc][r] for r in ResNames]
                    + [User_Angles[struc][f] for f in FacNames]
                )
                All_Angles[struc] = User_Angles[struc]
                All_Sequences[struc] = User_Sequences[struc]
            except KeyError, e:
                raise Exception("Structure %s not known" % e)

        try:
            combinations[comb] += 1
        except KeyError:
            combinations[comb] = 1

        for r in ResNames:
            rf[r][All_Sequences[struc][r]] += 1
        for f in FacNames:
            try:
                ff[f][All_Angles[struc][f]] += 1
            except KeyError:
                ff[f][All_Angles[struc][f]] = 1

    # Calculate the expected number of times we should see the combination based on the individual aa frequencies.
    expected = {}
    N = len(Names)
    nr = len(ResNames)
    nf = len(FacNames)
    n = nr + nf
    for comb in combinations:
        freq = [rf[ResNames[r]][comb[r]] for r in range(nr)] + [
            ff[FacNames[f]][comb[nr + f]] for f in range(nf)
        ]
        expected[comb] = int(round(float(reduce(mul, freq)) / (N ** (n - 1))))
    return combinations, expected, Set.Label, ResNames, FacNames


def CombinationSearch(args):
    """Find the possible combinations of amino acids and a given residue position(s)
    and output to stdout."""
    rsearch = ParseW(" ".join(args.r))[0]
    wasearch = ParseW(" ".join(args.wa))[0]
    fsearch = ParseW(" ".join(args.f))[0]
    wpsearch = ParseW(" ".join(args.wp))[0]
    out = sys.stdout
    if len(rsearch) <= len(wasearch) and len(fsearch) <= len(wpsearch):
        return False
    combinations, expected, label, rsearch, fsearch = GetFactorCombinations(
        rsearch, wasearch, fsearch, wpsearch, args
    )
    # Output to stdout. Header is not outputted if quiet option True.
    if not args.q:
        out.write(
            "Combinations and frequencies of requested residues and factors present"
        )
        if label == "all":
            out.write("\nin all structures\n")
        else:
            out.write("\nin structures with %s \n" % label)
        table = [rsearch + fsearch + ["Obs", "Expected"]]

    for c in combinations:
        table.append(list(c) + [combinations[c], expected[c]])
    dataIO.print_table(table, out)
    return True


def StructureSearch(args):
    """Take the query arguments and find structures which satisfy them"""

    # Parse the -r argument.
    ResidueArgs = ParseW(" ".join(args.r))

    # Parse the -wa ("with amino acids") argument.
    WithAAArgs = ParseW(args.wa)

    # The -r argument must either have only one value OR as many comma-separated values as the -wa argument
    nr = len(ResidueArgs)
    nwa = len(WithAAArgs)
    if nr == 1:
        ResidueArgs = ResidueArgs * nwa
    elif nr != nwa:
        raise Exception(
            'The number of comma separated sets of residues (-r argument) must be 1 OR as many comma-separated sets of "with amino-acid" -wa argument\n'
        )

    # Check that the length of each of the comma separated arguments in -r and -wa are the same.
    for i in range(nwa):
        if len(ResidueArgs[i]) != len(WithAAArgs[i]):
            raise Exception(
                "The number of amino acids specified is not the same as the number of residues specified in the set: %s ; %s\n"
                % (" ".join(ResidueArgs[i]), " ".join(WithAAArgs[i]))
            )

    # Parse the -f argument.
    FactorArgs = ParseW(" ".join(args.f))

    # Parse the -wp ("with parameter") argument.
    WithPArgs = ParseW(args.wp)

    # The -f argument must either have only one value OR as many comma-separated values as the -wp argument
    nf = len(FactorArgs)
    nwp = len(WithPArgs)
    if nf == 1:
        FactorArgs = FactorArgs * nwp
    elif nf != nwp:
        raise Exception(
            'The number of comma separated sets of factors (-f argument) must be 1 OR as many comma-separated sets of "with parameters" (-wp argument)\n'
        )

    # Check that the length of each of the comma separated arguments in -r and -wa are the same.
    for i in range(nwp):
        if len(FactorArgs[i]) != len(WithPArgs[i]):
            raise Exception(
                "The number of parameters specified is not the same as the number of factors specified in the set: %s ; %s\n"
                % (" ".join(FactorArgs[i]), " ".join(WithPArgs[i]))
            )

    if args.wa and args.wp:
        if nwp != nwa:
            if nwa == 1:
                ResidueArgs = ResidueArgs * nwp
                WithAAArgs = WithAAArgs * nwp
            elif nwp == 1:
                FactorArgs = FactorArgs * nwa
                WithPArgs = WithPArgs * nwa
            else:
                raise Exception(
                    "The number of comma-separated sets defined in the factor arguments is not the same as in the residue arguments\n"
                )

    # For each of the sets find the structures.
    SetsFound = []
    if args.wa and args.wp:
        for i in range(len(WithAAArgs)):
            SetsFound.append(
                StructureSet(ResidueArgs[i], WithAAArgs[i], FactorArgs[i], WithPArgs[i])
            )
    elif args.wa:
        for i in range(len(WithAAArgs)):
            SetsFound.append(StructureSet(ResidueArgs[i], WithAAArgs[i], [], []))
    elif args.wp:
        for i in range(len(WithPArgs)):
            SetsFound.append(StructureSet([], [], FactorArgs[i], WithPArgs[i]))

    return SetsFound


def UserSets(Usets):
    """Interpret the selection of user defined sets of structures from the -s argument"""
    Sets = []
    for s in ParseW(" ".join(Usets)):
        if s[0] == "all":
            s = Angles.keys()
            label = "pre-calculated structures"
        elif s[0] == "mine":
            if User_Angles:
                s = User_Angles.keys()
                label = "user's structures"
            else:
                sys.stderr.write("Warning: No stored user structures found\n")
                continue
        elif s[0].endswith(".dat"):
            # Expect a file  with a list of structure names. One per line.
            label = "structures defined in %s" % s[0]
            s = map(str.strip, open(s[0]).readlines())
        else:
            label = ""
        Sets.append(StructureSet([], [], [], [], s, label))
    return Sets


###########################
# Visualisation functions #
###########################


def OutputData(StructureSets, args):
    """Output the information for the structures requested"""
    out = sys.stdout
    for Set in StructureSets:
        table = []
        if not args.q:
            out.write("Structure(s) with:\n")
            out.write(" ".join(Set.Label.split("\n")) + "\n")
            table.append(["Structure"] + anglenames)

        for s in Set.structures:
            try:
                table.append([s] + [User_Angles[s][a] for a in anglenames])
            except KeyError:
                table.append([s] + [All_Angles[s][a] for a in anglenames])
        # Print list of lists to screen nicely
        dataIO.print_table(table, out)
        out.write("\n")


def ShowMsa(StructureSets, args):
    """Show the multiple sequence alignment of the structure sets selected over the residues requested"""
    if args.mr == "all":
        args.mr = Residues
    out = sys.stdout
    for Set in StructureSets:
        if not args.q:
            out.write("Structures with:\n" + Set.Label + "\n")
            table = [["Structure"] + args.mr]
        for s in Set.structures:

            try:
                table.append([s] + [User_Sequences[s][r] for r in args.mr])
            except KeyError:
                table.append([s] + [All_Sequences[s][r] for r in args.mr])
        # Print list of lists to screen nicely
        dataIO.print_table(table, out)
        out.write("\n")


def PlotAngles(StructureSets, args, NewAngles=[], show=True):
    """Plot the angles of the structures selected against the non-redundant background distributions.
    Plotting is done using an R script. Displaying the resultant plot is done using the program 'display' from ImageMagick suite.
    Comes with ubuntu linux distro."""

    Data = {}
    if args.plot == "all" or not args.plot:
        Data["anglenames"] = ["HL", "dc", "LC1", "HC1", "LC2", "HC2"]
    else:
        Data["anglenames"] = [args.plot]

    ######################################
    # Find the background distributions  #
    ######################################

    for ang in Data["anglenames"]:
        Data[ang + "Full"] = []
        for code in Angles.keys():
            Data[ang + "Full"].append("%.2f" % Angles[code][ang])

    ###############################
    # Allow for saving png output #
    ###############################
    if args.png:
        Data["fname"] = [args.png]
    else:
        Data["fname"] = [tempfile.mkstemp(".png", "Plot")[1]]

    ##############################
    # Set plot type for each set #
    ##############################
    Data["plottype"] = []
    Data["cols"] = []
    Data["Ltype"] = []
    colors = [
        "red",
        "blue",
        "green",
        "cyan",
        "orange",
        "purple",
        "pink",
    ]

    i = 0
    j = 1
    for sset in StructureSets:
        if not len(sset):
            continue
        if not args.lty:
            if len(sset) < 10:
                Data["plottype"].append(
                    "lines"
                )  # if there is only a few structures in the set
            else:  # default to density plot
                Data["plottype"].append("density")
        # else use what the user has asked for regardless of the set size.
        elif args.lty == "d":
            Data["plottype"].append("density")
        elif args.lty == "l":
            Data["plottype"].append("lines")
        sset.color = colors[i]
        Data["cols"].append(sset.color)
        Data["Ltype"].append("%d" % j)
        i += 1
        if i > len(colors) - 1:
            i = 0
            j += 1

    for n in NewAngles:
        NewAngles[n]["color"] = colors[i]
        Data["plottype"].append("lines")
        Data["cols"].append(colors[i])
        Data["Ltype"].append("2")
        i += 1
        if i > len(colors) - 1:
            i = 0

    #############################################################
    # Get the angle distribution for each of the structure sets #
    #############################################################

    for sset in StructureSets:
        if not len(sset):
            continue
        sset.GetAngles()

    n = 1
    for i in range(len(StructureSets)):
        if not len(StructureSets[i]):
            continue

        Data.update(
            dict(
                zip(
                    [a + str(n) for a in Data["anglenames"]],
                    [
                        map(lambda t: "%.2f" % t, StructureSets[i].Angles[a])
                        for a in Data["anglenames"]
                    ],
                )
            )
        )
        Data["Label" + str(n)] = (
            StructureSets[i]
            .Label.replace(" ", "@")
            .replace("\n", "/")
            .replace("'", "~")
        )
        n += 1
    for i in NewAngles:
        Data.update(
            dict(
                zip(
                    [a + str(n) for a in Data["anglenames"]],
                    [
                        map(lambda t: "%.2f" % t, [NewAngles[i][a]])
                        for a in Data["anglenames"]
                    ],
                )
            )
        )
        Data["Label" + str(n)] = (
            i.replace(" ", "@").replace("\n", "/").replace("'", "~")
        )
        n += 1

    ###############################
    # Put together the pipestring #
    ###############################
    pipestring = ""
    Fields = Data.keys()
    pipestring += "\t".join(Fields) + "\n"
    # Write to the pipestring as long as the longest vector nans otherwise.
    i = 0
    while 1:
        wrote = False
        for f in Fields:
            try:
                pipestring += Data[f][i] + "\t"
                wrote = True
            except IndexError:
                pipestring += "nan\t"
        pipestring += "\n"
        i += 1
        if not wrote:
            break

    if show:
        p = subprocess.Popen(
            ["Rscript", os.path.join(path, "plotting/plotting.R")],
            stdout=subprocess.PIPE,
            stdin=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        output = p.communicate(input=pipestring)
        if not args.png:
            os.system("display %s &" % Data["fname"][0])
            if StructureSets or NewAngles:
                os.system(
                    "display %s &" % Data["fname"][0].replace(".png", "legend.png")
                )

    if args.png:
        return []
    else:
        return [Data["fname"][0], Data["fname"][0].replace(".png", "legend.png")]


def ShowPymol(StructureSets, args):
    """Show structures in the sets objects listed in StructureSets. A maximum of 10 structures is shown from each set by default. Modify with command line option nalign"""
    import pymol
    from pymol import cmd, CmdException

    pymol.finish_launching()

    if args.nalign == "all":
        n = 10000
    else:
        try:
            n = int(args.nalign)
        except ValueError, e:
            sys.stderr.write("Warning: %s, setting nalign=10\n" % e)
            n = 10

    cmd.set(
        "bg_rgb", [1, 1, 1]
    )  # This sets the background to white - remove or change to [0,0,0] for black
    for Set in StructureSets:
        random.shuffle(Set.structures)
        for s in Set.structures[
            :n
        ]:  # Shows a maximum of 10 structures from each set by default
            name = s + "_" + Set.color
            try:
                if User_Angles.has_key(s):
                    cmd.load(
                        os.path.join(user_datapath, "user_fvs", s + ".pdb"), object=name
                    )
                else:
                    cmd.load(
                        os.path.join(data_path, "choth_fv", s + ".pdb"), object=name
                    )
            except CmdException:
                raise Exception("Structure File for %s not found\n" % s)
            cmd.hide("everything", name)
            align.tmalignconsensus(name, useH=args.useH)
            if args.showaxes:
                align.buildaxes(name, vectorlength=10, color=Set.color)

                cmd.set("dash_gap", 0.4)
                cmd.hide("labels")
            cmd.color(Set.color, name)

    cmd.show("lines", "  name c+ca+n")
    cmd.orient()
