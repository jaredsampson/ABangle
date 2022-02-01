from typing import List
import argparse
import pytest

from abangle.calculate import *
from abangle.number import *

@pytest.mark.parametrize(
        "fname,expected_angles", 
        [('4kq3.pdb', {'HC2':114.97, 'HC1':71.58, 'LC2':83.15, 'LC1':119.49, 'dc':16.00, 'HL':-61.10}), 
        ('2atk.pdb', {'HL': -54.09, 'HC1': 70.76, 'HC2': 114.09, 'LC1': 123.29, 'LC2': 83.42, 'dc': 16.48}),
        ('1u8l.pdb', {'HL': -56.41, 'HC1': 68.73, 'HC2': 118.87, 'LC1': 123.11, 'LC2': 82.99, 'dc': 15.88}
 )])
def test_find_angles(shared_datadir, fname, expected_angles):
    """Checks that angles computed by find_angles function are correct"""
    return {key:round(val, 2) for key, val in find_angles(shared_datadir/fname).items()}


if __name__ == '__main__':
    
    if False:
        examples = data_path.parent/'data'/'example_pdbs'

        angles_4KQ3 = {'HC2':114.97, 'HC1':71.58, 'LC2':83.15, 'LC1':119.49, 'dc':16.00, 'HL':-61.10}
        validate_angles('4kq3.pdb', **angles_4KQ3)

        angles_2ATK = {'HL': -54.09, 'HC1': 70.76, 'HC2': 114.09, 'LC1': 123.29, 'LC2': 83.42, 'dc': 16.48}
        validate_angles('2atk.pdb', **angles_2ATK)

        angles_1U8L = {'HL': -56.41, 'HC1': 68.73, 'HC2': 118.87, 'LC1': 123.11, 'LC2': 82.99, 'dc': 15.88}
        validate_angles('1u8l.pdb', **angles_1U8L)

