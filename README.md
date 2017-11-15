molecool
========
[![Build Status](https://travis-ci.org/steinmanngroup/molecool.svg?branch=master)](https://travis-ci.org/steinmanngroup/molecool)

Basic module to handle molecules which are represented as
a collection of atoms. A basic use-case example might be
to read all atoms from an .xyz-file using the Molecule
class to store atoms which are represented by the Atom
class.

```python
>>> with open('water.xyz', 'r') as f:
>>>     lines = f.readlines()
>>>
>>> molecule = Molecule()
>>> for atom_xyz in lines[2:]:
>>>     Z = molecule.util.LABEL2Z[atom_xyz[0]]
>>>     molecule.addAtom( Atom(Z, xyz=map(float, atom_xyz[1:])) )
>>>
>>> assert 3 == molecule.getNumAtoms()
```

The package also includes a wrapper around the openbabel
molecule class so it can be used as a backend. Hopefully
this will in the future be of some convenience instead
of me having to program SMARTS myself. Ugh.
