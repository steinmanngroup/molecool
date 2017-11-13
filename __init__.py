from atom import Atom
from molecule import Molecule

__author__ = "Casper Steinmann"
__copyright__ = "Copyright 2017"
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Casper Steinmann"
__email__ = "casper.steinmann@gmail.com"
__status__ = "Alpha"
__doc__ = """
Molecule
========

Basic module to handle molecules which are represented as
a collection of atoms. A basic use-case example might be
to read all atoms from an .xyz-file using the Molecule
class to store atoms which are represented by the Atom
class.

>>> with open('water.xyz', 'r') as f:
>>>     lines = f.readlines()
>>>
>>> molecule = Molecule()
>>> for atom_xyz in lines[2:]:
>>>     Z = molecule.util.LABEL2Z[atom_xyz[0]]
>>>     molecule.addAtom( Atom(Z, xyz=map(float, atom_xyz[1:])) )
>>>
>>> assert 3 == molecule.getNumAtoms()

"""
