'''Crystal definitions and SOAP vector calculations for simple
elements.
'''
# import numpy as np
from pyrelate import descriptors
from ase import Atoms


elements = {
    "Ni": ("FaceCenteredCubic", 3.52, 28, [0]),
    "Al": ("FaceCenteredCubic", 4.05, 13, [0]),
    "Cr": ("BodyCenteredCubic", 2.91, 24, [0, 1]),
    "Mg": ("HexagonalClosedPacked", {'a': 3.21, 'c/a': 1.633}, 12, [0, 1])
}
"""dict: keys are element names, values are a tuple of (`str` lattice,
`float` lattice parameter, `int` element number, `list` basis indices).
"""


def atoms(element):
    '''Returns an :class:`ase.Atoms` object for the given element, using the tabulated lattice parameters.

    Parameters:
        element (str): name of the element.
    '''
    lattice = "unknown"
    if element in elements:
        lattice, latpar, Z, basis = elements[element]
        if lattice == "HexagonalClosedPacked":
            import ase.lattice.hexagonal as structures
        else:
            import ase.lattice.cubic as structures
        if hasattr(structures, lattice):
            lat = getattr(structures, lattice)(element, latticeconstant=latpar)
            a = Atoms(positions=lat.positions, numbers=lat.numbers, pbc=True)
            a.set_cell(lat.cell)
            a.set_atomic_numbers([Z for i in a])
            return a
    else:
        raise NotImplementedError(f"Lattice info for element {element} not hard-coded into elements.py.")
    # FIXME throw error if structure not included


def seed(element, soap_fcn, **soapargs):
    """Computes the :math:`P` matrix for the given element.

    Parameters:
        element (str): name of the element.
        soap_fcn (function): function to compute SOAP matrix.
        soapargs (dict): Parameters associated with the SOAP description being used.
    """

    a = atoms(element)
    if soap_fcn is None:
        soap_fcn = descriptors.soap
    return soap_fcn(a, **soapargs)[0]
