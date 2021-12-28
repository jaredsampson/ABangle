from typing import List
from abangle.calculate import *
import argparse
import pytest
from abangle.calculate import *

@pytest.fixture
def heavy_sequence():
    return {
        'H2': 'V', 'H3': 'Q', 'H4': 'L', 'H5': 'V', 'H6': 'Q', 'H7': 'S', 'H8': 'G', 'H9': 'A', 'H10': 'E', 'H11': 'V', 'H12': 'K', 'H13': 'K', 
        'H14': 'P', 'H15': 'G', 'H16': 'S', 'H17': 'S', 'H18': 'V', 'H19': 'K', 'H20': 'V', 'H21': 'S', 'H22': 'C', 'H23': 'K', 'H24': 'A', 
        'H25': 'S', 'H26': 'G', 'H27': 'G', 'H28': 'T', 'H29': 'F', 'H30': 'S', 'H31': 'S', 'H32': 'Y', 'H33': 'A', 'H34': 'I', 'H35': 'S', 
        'H36': 'W', 'H37': 'V', 'H38': 'R', 'H39': 'Q', 'H40': 'A', 'H41': 'P', 'H42': 'G', 'H43': 'Q', 'H44': 'G', 'H45': 'L', 'H46': 'E', 
        'H47': 'W', 'H48': 'M', 'H49': 'G', 'H50': 'S', 'H51': 'I', 'H52': 'I', 'H52A': 'P', 'H53': 'W', 'H54': 'F', 'H55': 'G', 'H56': 'T', 
        'H57': 'T', 'H58': 'N', 'H59': 'Y', 'H60': 'A', 'H61': 'Q', 'H62': 'K', 'H63': 'F', 'H64': 'Q', 'H65': 'G', 'H66': 'R', 'H67': 'V', 
        'H68': 'T', 'H69': 'I', 'H70': 'T', 'H71': 'A', 'H72': 'D', 'H73': 'E', 'H74': 'S', 'H75': 'T', 'H76': 'S', 'H77': 'T', 'H78': 'A', 
        'H79': 'Y', 'H80': 'M', 'H81': 'E', 'H82': 'L', 'H82A': 'S', 'H82B': 'S', 'H82C': 'L', 'H83': 'R', 'H84': 'S', 'H85': 'E', 'H86': 'D', 
        'H87': 'T', 'H88': 'A', 'H89': 'V', 'H90': 'Y', 'H91': 'Y', 'H92': 'C', 'H93': 'A', 'H94': 'R', 'H95': 'D', 'H96': 'S', 'H97': 'E', 
        'H98': 'Y', 'H99': 'Y', 'H100': 'F', 'H101': 'D', 'H102': 'H', 'H103': 'W', 'H104': 'G', 'H105': 'Q', 'H106': 'G', 'H107': 'T', 
        'H108': 'L', 'H109': 'V', 'H110': 'T', 'H111': 'V', 'H112': 'S', 'H113': 'S'
    }

@pytest.fixture
def light_sequence():
    return {
        'L1': 'D', 'L2': 'I', 'L3': 'Q', 'L4': 'M', 'L5': 'T', 'L6': 'Q', 'L7': 'S', 'L8': 'P', 'L9': 'S', 'L10': 'S', 'L11': 'V', 'L12': 'S', 
        'L13': 'A', 'L14': 'S', 'L15': 'V', 'L16': 'G', 'L17': 'D', 'L18': 'R', 'L19': 'V', 'L20': 'T', 'L21': 'I', 'L22': 'T', 'L23': 'C', 
        'L24': 'R', 'L25': 'A', 'L26': 'S', 'L27': 'Q', 'L28': 'G', 'L29': 'I', 'L30': 'S', 'L31': 'N', 'L32': 'W', 'L33': 'L', 'L34': 'N', 
        'L35': 'W', 'L36': 'Y', 'L37': 'Q', 'L38': 'Q', 'L39': 'K', 'L40': 'P', 'L41': 'G', 'L42': 'K', 'L43': 'A', 'L44': 'P', 'L45': 'K', 
        'L46': 'L', 'L47': 'L', 'L48': 'I', 'L49': 'Y', 'L50': 'A', 'L51': 'A', 'L52': 'S', 'L53': 'S', 'L54': 'L', 'L55': 'Q', 'L56': 'S', 
        'L57': 'G', 'L58': 'V', 'L59': 'P', 'L60': 'S', 'L61': 'R', 'L62': 'F', 'L63': 'S', 'L64': 'G', 'L65': 'S', 'L66': 'G', 'L67': 'S', 
        'L68': 'G', 'L69': 'T', 'L70': 'D', 'L71': 'F', 'L72': 'T', 'L73': 'L', 'L74': 'T', 'L75': 'I', 'L76': 'S', 'L77': 'S', 'L78': 'L', 
        'L79': 'Q', 'L80': 'P', 'L81': 'E', 'L82': 'D', 'L83': 'F', 'L84': 'A', 'L85': 'T', 'L86': 'Y', 'L87': 'Y', 'L88': 'C', 'L89': 'Q', 
        'L90': 'Q', 'L91': 'Y', 'L92': 'S', 'L93': 'D', 'L94': 'D', 'L96': 'P', 'L97': 'T', 'L98': 'F', 'L99': 'G', 'L100': 'Q', 'L101': 'G', 
        'L102': 'T', 'L103': 'K', 'L104': 'V', 'L105': 'E', 'L106': 'I', 'L107': 'K', 'L108': 'R', 'L109': 'T', 'L110': 'V'}

def test_get_loop_length_H3(heavy_sequence):
    assert get_loop_length(heavy_sequence, 'H3') == 8

def test_get_loop_length_L1(light_sequence):
    assert get_loop_length(light_sequence, 'L1') == 11

def test_get_loop_length_with_gap(heavy_sequence):
    heavy_sequence.update({'H97': '-'})
    assert get_loop_length(heavy_sequence, 'H3') == 7

def test_get_loop_length_with_letters(heavy_sequence):
    heavy_sequence.update({'H97A': 'Y', 'H97B': 'H'})
    assert get_loop_length(heavy_sequence, 'H3') == 10