#!/usr/bin/env python3
import argparse
import csv
import json
import pathlib
import pprint
import sys

from abangle import calculate

description = """
DESCRIPTION	

	ABangle

	A tool to calculate and analyse the orientation between the VH and VL domains in antibodies.
	Command line interface.
	See README file for full description and examples.
"""

epilogue = """
AUTHOR
	
	2012
	James Dunbar - Oxford protein informatics group.
	james.dunbar@dtc.ox.ac.uk
	http://opig.stats.ox.ac.uk/
	http://www.stats.ox.ac.uk/~dunbar/
		
	Supervisors
	Prof Charlotte M Deane - Oxford protein informatics group.
	Dr Angelika Fuchs (Roche) and Dr Jiye Shi (UCB Celltech)
"""


def read_user_cli_args():
    """Handles the CLI user interactions."""
    parser = argparse.ArgumentParser(
        prog="ABangle",
        description=description,
        epilog=epilogue,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # For angle calculation for new structures from user
    parser.add_argument("-i", type=str, nargs="+", default=[], help="File name of structure to domain orientations for. Can also be a file containing a list of pdb files - one per line")
    parser.add_argument("-store", type=str, default="stdout", help="Specify format for saving angles data, by default result will be printed to stdout. User can specify data be stored in json or csv format")
    parser.add_argument("-outdir", type=str, default=".", help="Specify directory to store program output")
    parser.add_argument("-name", type=str, default="angle_data", help="The final path component of the output file, without its suffix. The suffix will be automatically determined from the '-store' argument")
    parser.add_argument("-scfv", type=str, nargs="+", default=[], help="If the molecule is in single chain fv format, use this argument to tell the program which chain identifier it is (or they are)")

    # For residue queries and finding structures
    parser.add_argument("-r", type=str, nargs="+", default=[], help="Specify residue positions. e.g. L44 H22 L2.")
    parser.add_argument("-wa", type=str, nargs="+", default=[], help="With amino acid. eg. A S D")
    parser.add_argument("-f", type=str, nargs="+", default=[], help="Factor query field. Choose from: HL, HC1, LC1, HC2, LC2, dc, Method, Res, Rfac, Bstate, AGtype, Hgroup, Lgroup, Ltype, Species, H3length.")
    parser.add_argument("-wp", type=str, nargs="+", default=[], help="Factor with parameters. For class factors put the class parameter e.g. -f Species -wp Mouse. For numeric factor use e.g: -f HL -wp dx:X   for HL within dx of x. Or -f HL -wp gt:X for HL greater than X")
    parser.add_argument("-s", type=str, nargs="+", default=[], help='Request individual structures. These can either be the PDB code e.g. 12E8  or the code and specific chain identifiers e.g. 12E8_HL. If it has been pre-calculated it will be reported. Use "mine" to define the users saved structures or "all" to define all pre-calculated structures')
    parser.add_argument("-target", type=str, default=False, help="Find Fvs which are structurally similar in terms of their angular measures to a target structure.")

    # Search set selection
    parser.add_argument("-useall", action="store_true", default=False, help="Use all the structures to search over. Not just a non-redundant set.")

    # Show all the info
    parser.add_argument("-showinfo", action="store_true", help="Show all the information for the structures in the set. Use -msa or -mr to see sequence")

    # For plotting measures
    parser.add_argument("-plot", type=str, default=False, choices=["all", "HL", "HC1", "HC2", "LC1", "LC2", "dc"], help="Show a plot of the angle distributions of selected sets of structures against background distribution. Use 'all' to plot all angles")
    parser.add_argument("-lty", type=str, choices=["l", "d"], default="", help="Show the sets using a density plot or vertical lines. Default makes a decision based on set size")
    parser.add_argument("-bg", type=str, default=False, help="File containing the data for backgound required for plots. Default is a non-redundant set.")
    parser.add_argument("-png", type=str, default=False, help="Save the angle plots as .png image file. No output to screen")

    # For showing structures in pymol
    parser.add_argument("-pymol", action="store_true", help="Launch pymol and overlay selected structures. Structures are all aligned to the L chain core consensus strucuture.")
    parser.add_argument("-nalign", type=str, default="10", help='The maximum number of structures from each set to visualise in pymol. Use "all" to show all (Use with caution for large sets)')
    parser.add_argument("-useH", action="store_true", default=False, help="Align the structures to the VH consensus domain in pymol")
    parser.add_argument("-showaxes", action="store_true", default=False, help="Show the axes in pymol")

    # For sequence or other information output
    parser.add_argument("-msa", action="store_true", help="Output msa over Fv region")
    parser.add_argument("-mr", type=str, nargs="+", default=False, help="Show msa of selected residues")
    parser.add_argument("-seqid", action="store_true", help="Output the top 10 highest matched sequence identity pairs of Fvs from the first two sets defined (or within the first set if only one)")

    return parser.parse_args()


def is_pdb_file(f):
    return pathlib.Path(f).suffix == ".pdb"


def as_posixpaths(files):
    return [pathlib.Path(f) for f in files]


def inputs_listed_in_file(user_args):
    return len(user_args.i) == 1 and not is_pdb_file(user_args.i[0])


def read_input_files(user_args):
    f = pathlib.Path(user_args.i[0])
    input_files = f.read_text().splitlines()
    assert all(is_pdb_file(f) for f in input_files), "some files may not be valid structure files"

    return input_files


def write_angles_to_csv(angles, outdir, name):
    with open((outdir / name).with_suffix(".csv"), "w", newline="") as f:
        fieldnames = ["pdb_id", "HL", "LC1", "LC2", "HC1", "HC2", "dc"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for pdb_id, angle in angles.items():
            row = dict(angle)
            row["pdb_id"] = pdb_id
            writer.writerow(row)


def write_angles_to_json(angles, outdir, name):
    with open((outdir / name).with_suffix(".json"), "w") as f:
        json.dump(angles, f)


def main():
    user_args = read_user_cli_args()

    if inputs_listed_in_file(user_args):
        input_files = read_input_files(user_args)
    else:
        input_files = user_args.i

    paths = as_posixpaths(input_files)
    angles = {path.stem: calculate.find_angles(path) for path in paths}

    outdir = pathlib.Path(user_args.outdir)
    if not outdir.exists():
        outdir.mkdir(parents=True)
    name = pathlib.Path(user_args.name)

    if user_args.store == "stdout":
        pprint.pprint(angles)
    elif user_args.store == "json":
        write_angles_to_json(angles, outdir, name)
    elif user_args.store == "csv":
        write_angles_to_csv(angles, outdir, name)
    else:
        raise ValueError("Storage format not supported")


if __name__ == "__main__":
    sys.exit(main())
