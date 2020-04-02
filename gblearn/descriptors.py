from pycsoap.soaplite import SOAP
import numpy as np

def soap(atoms, rcut, nmax, lmax, **kwargs):
    """
        nmax (int): bandwidth limits for the SOAP descriptor radial basis
          functions.
        lmax (int): bandwidth limits for the SOAP descriptor spherical
          harmonics.
        rcut (float): local environment finite cutoff parameter.
    """
    soap_desc = SOAP(atomic_numbers=atoms.numbers, rcut=rcut,
                     nmax=nmax, lmax=lmax, **kwargs)
    P = soap_desc.create(atoms)
    return P

def asr(atoms, store, rcut, nmax, lmax, **kwargs):
    aid = atoms.get_array("aid")[0]
    matrix = store.get(
        "soap", aid, rcut=rcut, nmax=nmax, lmax=nmax, **kwargs)
    #the collection is written into the filename, we don't really want that...
    return np.average(matrix, axis=0)

def ler(atoms, collection, eps, rcut, nmax, lmax, seed=None, **kwargs):
    """
    """
    U = collection.store.get("ler", 'U', rcut=rcut, nmax=nmax, lmax=lmax, **kwargs)
    if U is None:
        from collections import OrderedDict
        U = {
            'centers': OrderedDict(),
            'clusters': {}
        }

        #Part 1: Clustering
        import gblearn.elements as elements
        #i added the extra args so that the FCC has the same params, was that right?
        U['centers'][('0', 0)] = elements.seed(collection[0].get_atomic_symbols()[0], rcut=rcut, nmax=nmax, lmax=lmax, **kwargs)
        for aid in collection:
            for lae_num, lae in enumerate(collection.store.get("soap", aid, rcut=rcut, nmax=nmax, lmax=lmax, **kwargs)):
                if lae is None:
                    raise RuntimeError("LER requries SOAP to be generated first")
                for unique in U['centers']:
                    dist = np.linalg.norm(unique - lae)
                    if dist < eps:
                        break
                else:
                    U['clusters'][(aid,lae_num)] = np.copy(lae)

        #Part 2: Classifying
        import annoy


def test(self, **args):
    return 'test result'
