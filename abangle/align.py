"""
DESCRIPTION	
	
	abangle.align.py
	A module to align structures in PyMOL using tmalign. 

OUTPUT
	
	~ Function to aligning Fv regions to the VL or VH consensus structures
	~ Function to show the ABangle coordinate system registered onto each Fv structure

AUTHOR
	
	2012
	James Dunbar - Oxford protein informatics group.
	james.dunbar@dtc.ox.ac.uk
	www.stats.ox.ac.uk/~dunbar/
		
	Supervisors
	Prof C.Deane - Oxford protein informatics group.
	Dr Angelika Fuchs (Roche) and Dr Jiye Shi (UCB Celltech)\n
"""

try:
    from pymol import cmd, CmdException
except ImportError:
    pass
import subprocess, tempfile, os, re
from pathlib import Path
from abangle import calculate

data_path = str(Path(__file__).parents[1]/'data')
coresetL = [
    int(l[1:])
    for l in map(str.strip, open(os.path.join(data_path, "Lcoresetfw.txt")).readlines())
]
coresetH = [
    int(l[1:])
    for l in map(str.strip, open(os.path.join(data_path, "Hcoresetfw.txt")).readlines())
]


def save_pdb_without_ter(filename, selection, **kwargs):
    """
    DESCRIPTION
        Save PDB file without TER records
    """
    v = cmd.get_setting_boolean("pdb_use_ter_records")
    if v:
        cmd.unset("pdb_use_ter_records")
    cmd.save(filename, selection, **kwargs)
    if v:
        cmd.set("pdb_use_ter_records")


def tmalignconsensus(
    mobile, args="", exe="TMalign", transform=1, object=None, quiet=1, useH=False
):
    quiet = int(quiet)

    mobile_filename = tempfile.mktemp(".pdb", "mobile")
    if useH:
        target_filename = os.path.join(data_path, "consensus_H.pdb")
        mobile_ca_sele = (
            "(%s) and (not hetatm) and name CA and alt +A and chain H and resi "
            % (mobile)
        )
        for res in coresetH:
            mobile_ca_sele += "%d+" % res
    else:
        target_filename = os.path.join(data_path, "consensus_L.pdb")
        mobile_ca_sele = (
            "(%s) and (not hetatm) and name CA and alt +A and chain L and resi "
            % (mobile)
        )
        for res in coresetL:
            mobile_ca_sele += "%d+" % res

    mobile_ca_sele = mobile_ca_sele[:-1]

    save_pdb_without_ter(mobile_filename, mobile_ca_sele)

    exe = cmd.exp_path(exe)
    matrix_filename = tempfile.mktemp(".txt", "matrix")
    args = [exe, mobile_filename, target_filename, "-m", matrix_filename] + args.split()

    try:
        process = subprocess.Popen(args, stdout=subprocess.PIPE)
    except OSError as exe:
        print(f'Cannot execute {exe}, please provide full path to TMalign executable')
        raise CmdException

    r = None
    re_score = re.compile(r"TM-score\s*=\s*(\d*\.\d*)")
    rowcount = 0
    matrix = []
    line_it = iter(process.stdout)
    attempt = 0
    while 1:
        for line in line_it:
            if 4 >= rowcount > 0:
                if rowcount >= 2:
                    a = map(float, line.split())
                    matrix.extend(a[2:5])
                    matrix.append(a[1])
                rowcount += 1
            elif line.lower().startswith(" -------- rotation matrix"):
                rowcount = 1
            elif line.startswith('(":" denotes'):
                alignment = [line_it.next().rstrip() for i in range(3)]
            else:
                match = re_score.search(line)
                if match is not None:
                    r = float(match.group(1))
            if not quiet:
                print(line.rstrip())
        if matrix or attempt:
            break
        else:
            try:
                line_it = iter(open(matrix_filename))
                attempt = 1
            except IOError:
                break
    if not quiet:
        for i in range(0, len(alignment[0]) - 1, 78):
            for line in alignment:
                print(line[i : i + 78])
            print('')

    assert len(matrix) == 3 * 4
    matrix.extend([0, 0, 0, 1])

    if int(transform):
        cmd.transform_selection("byobject (%s)" % (mobile), matrix, homogenous=1)

    # alignment object
    if object is not None:
        mobile_idx, target_idx = [], []
        space = {"mobile_idx": mobile_idx, "target_idx": target_idx}
        cmd.iterate(
            mobile_ca_sele, 'mobile_idx.append("%s`%d" % (model, index))', space=space
        )
        cmd.iterate(
            target_ca_sele, 'target_idx.append("%s`%d" % (model, index))', space=space
        )
        for i, aa in enumerate(alignment[0]):
            if aa == "-":
                mobile_idx.insert(i, None)
        for i, aa in enumerate(alignment[2]):
            if aa == "-":
                target_idx.insert(i, None)
        if len(mobile_idx) == len(target_idx) == len(alignment[2]):
            cmd.rms_cur(
                " ".join(
                    idx for (idx, m) in zip(mobile_idx, alignment[1]) if m in ":."
                ),
                " ".join(
                    idx for (idx, m) in zip(target_idx, alignment[1]) if m in ":."
                ),
                cycles=0,
                matchmaker=4,
                object=object,
            )
        else:
            print("Could not load alignment object")

    if not quiet and r is not None:
        print(f"Found in output TM-score = {r:.4f}")

    os.remove(mobile_filename)
    try:
        os.remove(matrix_filename)
    except OSError:
        pass
    return r


def buildaxes(name, vectorlength=10, color="yellow"):
    """Function to show the coordinate system in pymol registered onto a structure"""
    objectname = name
    vectorlength = float(vectorlength)
    _tmp_fname = tempfile.mkstemp(".pdb", "name")[1]
    cmd.save(_tmp_fname, objectname)

    points = {}
    Lpoints, Hpoints = calculate.mapvectors(_tmp_fname, [])

    C = normalise([Hpoints[0][i] - Lpoints[0][i] for i in range(3)])
    Cminus = map(lambda x: -1 * x, C)
    L1 = normalise([Lpoints[1][i] - Lpoints[0][i] for i in range(3)])
    L2 = normalise([Lpoints[2][i] - Lpoints[0][i] for i in range(3)])
    H1 = normalise([Hpoints[1][i] - Hpoints[0][i] for i in range(3)])
    H2 = normalise([Hpoints[2][i] - Hpoints[0][i] for i in range(3)])
    dc = (
        sum(map(lambda x: x ** 2, [Hpoints[0][i] - Lpoints[0][i] for i in range(3)]))
        ** 0.5
    )
    points["VLC"] = Lpoints[0]
    points["VHC"] = Hpoints[0]
    points["VL1"] = ADD(
        points["VLC"],
        normalise([Lpoints[1][i] - Lpoints[0][i] for i in range(3)], vectorlength),
    )
    points["VL2"] = ADD(
        points["VLC"],
        normalise([Lpoints[2][i] - Lpoints[0][i] for i in range(3)], vectorlength),
    )
    points["VH1"] = ADD(
        points["VHC"],
        normalise([Hpoints[1][i] - Hpoints[0][i] for i in range(3)], vectorlength),
    )
    points["VH2"] = ADD(
        points["VHC"],
        normalise([Hpoints[2][i] - Hpoints[0][i] for i in range(3)], vectorlength),
    )

    numbers = dict(
        zip(
            ["VLC", "VHC", "VL1", "VL2", "VH1", "VH2"],
            [1001, 1002, 1003, 1004, 1005, 1006],
        )
    )

    for p in points:
        cmd.pseudoatom(
            objectname,
            resi=numbers[p],
            pos=points[p],
            name="ca",
            color="green",
            hetatm=True,
        )

    cmd.color("red", "%s and resi %d+%d" % (objectname, numbers["VH1"], numbers["VL1"]))
    cmd.color(
        "pink", "%s and resi %d+%d" % (objectname, numbers["VLC"], numbers["VHC"])
    )

    cmd.distance(
        objectname + "dc",
        "%s and resi %d" % (objectname, numbers["VLC"]),
        "%s and resi %d" % (objectname, numbers["VHC"]),
    )
    cmd.distance(
        objectname + "L1",
        "%s and resi %d" % (objectname, numbers["VLC"]),
        "%s and resi %d" % (objectname, numbers["VL1"]),
    )
    cmd.distance(
        objectname + "L2",
        "%s and resi %d" % (objectname, numbers["VLC"]),
        "%s and resi %d" % (objectname, numbers["VL2"]),
    )
    cmd.distance(
        objectname + "H1",
        "%s and resi %d" % (objectname, numbers["VHC"]),
        "%s and resi %d" % (objectname, numbers["VH1"]),
    )
    cmd.distance(
        objectname + "H2",
        "%s and resi %d" % (objectname, numbers["VHC"]),
        "%s and resi %d" % (objectname, numbers["VH2"]),
    )

    cmd.color(color, objectname + "dc")
    cmd.color(color, objectname + "L1")
    cmd.color(color, objectname + "L2")
    cmd.color(color, objectname + "H1")
    cmd.color(color, objectname + "H2")

    cmd.show("nb_spheres", objectname)

    os.remove(_tmp_fname)


def normalise(a, length=1):
    mag = (sum(map(lambda x: x ** 2, a))) ** 0.5
    return map(lambda x: length * x / mag, a)


def ADD(a, b):
    return [a[i] + b[i] for i in range(3)]

if __name__ == '__main__':
    print(os.path.join(os.path.split(__file__)[0], "data"))