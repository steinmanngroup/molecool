import numpy
""" Utility variables and functions
"""

aa2au = 1.8897261249935897  # bohr / AA

# converts nuclear charge to atom label
Z2LABEL = {
 1: 'H',                                                             2: 'He',
 3: 'Li',  4: 'Be',  5: 'B',   6: 'C',   7: 'N',  8: 'O',  9: 'F',  10: 'Ne',
11: 'NA', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar'
}

# converts an atomic label to a nuclear charge
LABEL2Z = {}
for key in Z2LABEL:
    LABEL2Z[Z2LABEL[key]] = key

# masses from UIPAC: http://www.chem.qmul.ac.uk/iupac/AtWt/
MASSES = {
  0: 0.00,
  1:  1.00784,                                                                                2: 4.002602,
  3:  6.938,   4:  9.01218, 5: 10.806,   6: 12.0096, 7: 14.00643, 8: 15.99903, 9: 18.998403, 10: 20.1797,
 11: 22.9898, 12: 24.304,  13: 26.9815, 14: 28.084, 15: 30.973,  16: 32.059,  17: 35.446,    18: 39.948
}

# Van der Waal radii from Alvarez (2013), DOI: 2013/dt/c3dt50599e
# all values in Angstrom
VDWRADII = {0: 0.00,
  1: 1.20,                                                              2: 1.43,
  3: 2.12,  4: 1.98,  5: 1.91,  6: 1.77,  7: 1.66,  8: 1.50,  9: 1.46, 10: 1.58,
 11: 2.50, 12: 2.51, 13: 2.25, 14: 2.19, 15: 1.90, 16: 1.89, 17: 1.82, 18: 1.83
}

# Covalent radii from Pykko and Atsumi (2009), DOI: 0.1002/chem.200800987
# all values in Angstrom
COVALENTRADII = {0: 0.00,
  1: 0.32,                                                              2: 0.46,
  3: 1.33,  4: 1.02,  5: 0.85,  6: 0.75,  7: 0.71,  8: 0.63,  9: 0.64, 10: 0.67,
 11: 1.55, 12: 1.39, 13: 1.26, 14: 1.16, 15: 1.11, 16: 1.03, 17: 0.99, 18: 0.96
}

# Coordination numbers from Pykko and Atsumi (2009), DOI: 0.1002/chem.200800987
COORDINATION = {0: 0,
  1: 1,                                            2: 1,
  3: 1,  4: 2,  5: 3,  6: 4,  7: 3,  8: 2,  9: 1, 10: 1,
 11: 1, 12: 2, 13: 3, 14: 4, 15: 3, 16: 2, 17: 1, 18: 1
}


def idamax(a):
    """ Returns the index of maximum absolute value (positive or negative)
        in the input array a.

        Note: Loosely based of a subroutine in GAMESS with the same name

        Arguments:
        a -- a numpy array where we are to find the maximum
             value in (either positive or negative)

        Returns:
        the index in the array where the maximum value is.
    """
    idx = -1
    v = 0.0
    for i, value in enumerate(numpy.abs(a)):
        if value > v:
            idx = i
            v = value

    return idx

def idamin(a):
    """ Returns the index of minimum absolute value (positive or negative)
        in the input array a.

        Arguments:
        a -- a numpy array where we are to find the minimum
             value in (either positive or negative)

        Returns:
        the index in the array where the maximum value is.
    """
    idx = -1
    v = 1.0e30
    for i, value in enumerate(numpy.abs(a)):
        if value < v:
            idx = i
            v = value

    return idx

def file_to_obmol(filename):
    file_format = obformat_from_filename(filename)
    mol = obmol_from_filename_and_format(filename, file_format)
    return mol

def obformat_from_filename(filename):
    return file_extension(filename)[1:]

def obmol_from_filename_and_format(filename, file_format):
    obc = openbabel.OBConversion()
    obc.SetInFormat(file_format)
    mol = openbabel.OBMol()
    obc.ReadFile(mol, filename)
    return mol

def file_extension(path_to_file):
    (filename, extension) = get_filename_and_extension(path_to_file)
    return extension

def file_basename(path_to_file):
    (filename, extension) = get_filename_and_extension(path_to_file)
    return filename

def get_filename_and_extension(path_to_file):
    if not isinstance(path_to_file, str):
        raise TypeError
    basename = os.path.split(path_to_file)[1]
    return os.path.splitext(basename)
