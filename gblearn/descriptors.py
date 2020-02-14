from dscribe.descriptors import SOAP


def soap(atoms, species, rcut,nmax,lmax):
    #TODO add args for additional parameters
    rbf = "gto"
    soap_desc = SOAP(species=species, periodic=False, rcut=rcut, nmax=nmax, lmax=lmax)
    P = soap_desc.create(atoms)
    return P

def test(idd, species, rcut, nmax, lmax):
    import numpy as np
    arb = np.array([[1,1,1],[2,2,2]])
    return arb
