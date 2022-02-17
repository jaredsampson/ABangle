import json
import re

with open('data/Sequences.dat') as f:
    data = [line.split() for line in f]

# define framework ranges
frameworks_l = set(
    list(range(1,25)) + \
    list(range(33,50)) + \
    list(range(62,91)) + \
    list(range(97,108))
)

frameworks_h = set(
    list(range(1,24)) + \
    list(range(34,51)) + \
    list(range(67,95)) + \
    list(range(102,114))
)

def resnum(s):
    """Get the residue number from position ID e.g. H12 -> 12"""
    match = re.search('\d+', s)
    return int(match.group()) if match else match

def chain(s):
    """Get the chain from the position ID e.g. 'H12' -> 'H'"""
    match = re.search('^[HL]', s)
    return match.group() if match else (match) 

def is_framework(posid):
    """Checks if a position ID falls in framework region"""
    return (chain(posid) == 'H' and resnum(posid) in frameworks_h) or \
           (chain(posid) == 'L' and resnum(posid) in frameworks_l)

def reformat_seq(res_list):
    """Joins list of characters and removes gaps e.g. ['T', 'V', '-', 'Y'] -> 'TVY'"""
    return ''.join(res_list).replace('-', '')

codes = data[0] # position codes e.g H6, H6A, H6B
seqs = data[1:] # residues e.g. ['E', 'V', 'Q', ...]

# filters residues that lie in the framework region and concatenates them.
# result is dictionary mapping pdb id to dictionary of heavy and light framework sequences
heavy_fws, light_fws = {}, {}
for seq in seqs:
    res_list = []
    for posid, res in zip(codes, seq):
        if is_framework(posid):
            res_list.append(res)
    # elements :85 are heavy chain positions, elements 85: are light chain positions
    H, L = reformat_seq(res_list[:85]), reformat_seq(res_list[85:])
    heavy_fws[seq[0]] = H
    light_fws[seq[0]] = L
    
# dump results into json file
with open('data/heavy_fws.json', 'w') as f:
    json.dump(heavy_fws, f)

with open('data/light_fws.json', 'w') as f:
    json.dump(light_fws, f)
