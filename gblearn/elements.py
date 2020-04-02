"""Crystal definitions and SOAP vector calculations for simple
elements.
"""
import numpy as np
from gblearn import descriptors
from ase import Atoms
_shells = {}
"""dict: keys are element names, values are lists of shells (in Ang.).
"""

elements = {
    "Ni": ("FaceCenteredCubic", 3.52, 28, [0]),
    "Al": ("FaceCenteredCubic", 4.05, 13, [0]),
    "Cr": ("BodyCenteredCubic", 2.91, 24, [0, 1]),
    "Mg": ("HexagonalClosedPacked", {'a':3.21, 'c/a':1.633}, 12, [0, 1])
}
"""dict: keys are element names, values are a tuple of (`str` lattice,
`float` lattice parameter, `int` element number, `list` basis indices).
"""

def atoms(element):
    """Returns a :class:`quippy.Atoms` structure for the given
    element, using the tabulated lattice parameters.
    Args:
        element (str): name of the element.
    """
    lattice = "unknown"
    if element in elements:
        lattice, latpar, Z, basis = elements[element]
        if lattice == "HexagonalClosedPacked":
            import ase.lattice.hexagonal as structures
        else:
            import ase.lattice.cubic as structures
        if hasattr(structures, lattice):
            lat = getattr(structures, lattice)(element, latticeconstant=latpar)
            a = Atoms(positions=lat.positions, numbers=lat.numbers)
            a.set_cell(lat.cell)
            a.set_atomic_numbers([Z for i in a])
            return a
#do rcut and nmax and lmax need to match what we call in SOAP?
def seed(element, lmax=12, nmax=12, rcut=6.0, *kwargs):
    """Computes the :math:`P` matrix for the given element.
    Args:
        element (str): name of the element.
        nmax (int): bandwidth limits for the SOAP descriptor radial basis
          functions.
        lmax (int): bandwidth limits for the SOAP descriptor spherical
          harmonics.
        rcut (float): local environment finite cutoff parameter.
    """
    lattice, latpar, Z, basis = elements[element]
    a = atoms(element)
    return descriptors.soap(a, rcut=rcut, nmax=nmax, lmax=lmax, **kwargs)
