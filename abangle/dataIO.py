"""
DESCRIPTION	
	
	abangle.dataIO.py
	A module to load data in abangle

OUTPUT
	
	~ Data structure for abangle

AUTHOR
	
	2012
	James Dunbar - Oxford protein informatics group.
	james.dunbar@dtc.ox.ac.uk
	www.stats.ox.ac.uk/~dunbar/
		
	Supervisors
	Prof C.Deane - Oxford protein informatics group.
	Dr Angelika Fuchs (Roche) and Dr Jiye Shi (UCB Celltech)\n
"""

import sys


def load(
    fname,
    header=False,
    rownames=False,
    colnames=False,
    coltype=dict,
    rowtype=dict,
    conv=True,
    transpose=False,
):
    __help__ = """Function to load a data file.
	fname -  the file name and path.
	header - specify whether the column names are in a header
	rownames  - False if no row names. Either specify the column with which to use as row name or provide a string list
	rowtype - store the entries for the row in a dictionary or a list.
	coltype - store the entries for the column variables in a dictionary, list or string (transpose=False).
	conv - convert the variables to floats if possible.
	transpose - changes so that the returned data structure indexes by column to row instead. e.g dataout =  { 'var1': { 'row1':vr11} } 
	e.g. rowtype = dict, coltype = dict: dataout = { 'row1': { 'var1': vr11 , 'var2: vr12}, 'row2': { 'var1': vr21 , 'var2: vr22} }
	     rowtype = list, coltype = dict: dataout = [ { 'var1': vr11 , 'var2: vr12},  { 'var1': vr21 , 'var2: vr22} ]
		
	Values:
		if header: returns dataout and colnames
		else: return dataout
	"""

    datafile = open(fname).readlines()

    if conv:
        dat = [list(map(convert, l.split())) for l in datafile]
    else:
        dat = list(map(str.split, datafile))

    if header:
        colnames = dat[0]
        j = 1
    elif (type(colnames) is list) and (
        (len(colnames) == len(dat[1]))
        or ((len(colnames) == len(dat[1]) - 1) and type(rownames) is int)
    ):
        # colnames correct length
        j = 0
    elif colnames:
        sys.stderr.write("Error: colnames is not the correct length\n")
        raise Exception
    else:
        j = 0

    if type(rownames) is int:
        rownames = [d[rownames] for d in dat[j:]]
        i = 1
    elif (type(rownames) is list) and (len(rownames) == len(dat) - j):
        # rownames correct length
        i = 0
    elif rownames:
        sys.stderr.write(
            "Erroprint_tabler: rownames must be column number or list of length data\n"
        )
        raise Exception
    else:
        i = 0

    colnames = colnames[i:]
    if transpose:
        if coltype is dict:
            if rowtype is dict:
                dataout = dict(
                    list(zip(
                        colnames,
                        (
                            dict(
                                (rownames[u], dat[j:][u][i:][v])
                                for u in range(len(rownames))
                            )
                            for v in range(len(colnames))
                        ),
                    ))
                )
            if rowtype is list:
                dataout = dict(
                    list(zip(
                        colnames,
                        (
                            [dat[j:][u][i:][v] for u in range(len(rownames))]
                            for v in range(len(colnames))
                        ),
                    ))
                )
        else:  # coltype is list
            if rowtype is dict:
                dataout = [
                    dict((rownames[u], dat[j:][u][i:][v]) for u in range(len(rownames)))
                    for v in range(len(colnames))
                ]
            if rowtype is list:
                dataout = [
                    [dat[j:][u][i:][v] for u in range(len(rownames))]
                    for v in range(len(colnames))
                ]
    else:
        if rowtype is dict:
            if coltype is dict:
                dataout = dict(
                    list(zip(rownames, [dict(list(zip(colnames, x[i:]))) for x in dat[j:]]))
                )
            elif coltype is list:
                dataout = dict(list(zip(rownames, [x[i:] for x in dat[j:]])))
            elif coltype is str:
                dataout = dict(
                    (rownames[r], datafile[j:][r][len(rownames[r]) :])
                    for r in range(len(rownames))
                )
        elif rowtype is list:
            if coltype is dict:
                dataout = [dict(list(zip(colnames, x[i:]))) for x in dat[j:]]
            elif coltype is list:
                dataout = [x[i:] for x in dat[j:]]
            elif coltype is str:
                dataout = datafile[j:]

    if header:
        return dataout, colnames
    else:
        return dataout


def convert(x):
    try:
        x = float(x)
    except ValueError:
        pass
    return x


def print_table(table, out):
    """Format the table (list of lists) and print to file object out"""
    col_paddings = []
    for i in range(len(table[0])):
        col_paddings.append(get_max_width(table, i))
    for row in table:
        print(row[0].ljust(col_paddings[0] + 1), end=' ', file=out)
        for i in range(1, len(row)):
            col = format_data(row[i]).rjust(col_paddings[i] + 1)
            print(col, end=' ', file=out)
        print(file=out)


def format_data(dat):
    if type(dat) is float:
        return "%.2f" % dat
    elif type(dat) is str:
        return dat
    elif type(dat) is int:
        return "%d" % dat
    else:
        raise Exception("Formatting error")


def get_max_width(table, index):
    return max([len(format_data(row[index])) for row in table])

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
        dat["H3length"] = get_loop_length(seq, 'H3')
        dat["L1length"] = get_loop_length(seq, 'L1')
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
            seq.append(aa_code_dict[aa])  # try to append the one letter code to the sequence.
        except KeyError:  # if we fail then there may be multiple occupancy
            if len(aa) == 4:  # check for residue name such as ASER or BSER
                aa3 = aa[1:]
                try:
                    seq.append(aa_code_dict[aa3])
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
