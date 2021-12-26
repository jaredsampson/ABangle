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
