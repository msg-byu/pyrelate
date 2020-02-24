from pycsoap.soaplite import SOAP
def soap(atoms, atomic_numbers, rcut,nmax,lmax, **kwargs):
    #TODO add args for additional parameters
    soap_desc = SOAP(atomic_numbers=atomic_numbers, rcut=rcut, nmax=nmax, lmax=lmax)
    P = soap_desc.create(atoms)
    return P

def test(aid, species, rcut, nmax, lmax):
    import numpy as np
    arb = np.array([[1,1,1],[2,2,2]])
    return arb
