from dscribe.descriptors import SOAP
import collection as C
def soap(idd, species, rcut,nmax,lmax):#add args for additional parameters
    rbf = "gto"
    soap_desc = SOAP(species=species, periodic=False, rcut=rcut, nmax=nmax, lmax=lmax)
    P = soap_desc.create(C.Collection.atoms_files[idd])
    return P

def test(idd, species, rcut, nmax, lmax):
    import numpy as np
    arb = np.array([[1,1,1],[2,2,2]])
    return arb
