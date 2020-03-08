from pycsoap.soaplite import SOAP


def soap(atoms, atomic_numbers, rcut, nmax, lmax, **kwargs):
    soap_desc = SOAP(atomic_numbers=atomic_numbers, rcut=rcut,
                     nmax=nmax, lmax=lmax, **kwargs)
    P = soap_desc.create(atoms)
    return P


def test(self, **args):
    return 'test result'
