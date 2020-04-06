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
    return np.average(matrix, axis=0)

def ler(atoms, collection, eps, rcut, nmax, lmax, seed=None, metric='euclidean', n_trees=10, search_k=-1, **kwargs):
    U = collection.store.get(
        "ler", 'U', collection=collection, eps=eps, rcut=rcut, nmax=nmax, lmax=lmax, metric=metric, n_trees=n_trees, search_k=search_k, **kwargs)#add seed? or hash?
    if U is None:
        from collections import OrderedDict
        U = {
            'centers': OrderedDict(),
            'clusters': {}
        }

        #Part 1: Clustering
        import gblearn.elements as elements
        if seed is None:
            seed = elements.seed(list(collection.values())[0].get_chemical_symbols()[0], rcut=rcut, nmax=nmax, lmax=lmax, **kwargs)
        U['centers'][('0', 0)] = seed
        for aid in collection:
            for lae_num, lae in enumerate(collection.store.get("soap", aid, rcut=rcut, nmax=nmax, lmax=lmax, **kwargs)):
                if lae is None:
                    raise RuntimeError("LER requries SOAP to be generated first")
                for unique in U['centers'].values():
                    dist = np.linalg.norm(unique - lae)
                    if dist < eps:
                        break
                else:
                    U['centers'][(aid,lae_num)]= np.copy(lae)
        #FIXME sort the unique bins

        #Part 2: Classifying
        from annoy import AnnoyIndex
        dim = len(seed)
        a = AnnoyIndex(dim, metric)
        for i, unique in enumerate(U['centers'].values()):
            a.add_item(i, unique)
        a.build(n_trees)
        a.save('tmp') #TODO find better way to save temporary file
        for aid in collection:
            for lae_num, lae in enumerate(collection.store.get("soap", aid, rcut=rcut, nmax=nmax, lmax=lmax, **kwargs)):
                nni = a.get_nns_by_vector(lae, 1, search_k=search_k, include_distances=False)[0]
                cluster = list(U['centers'].keys())[nni]
                #FIXME when we sort just fill the empty lists first
                if cluster not in U['clusters']:
                    U['clusters'][cluster] = []
                U['clusters'][cluster].append((aid, lae_num))
        import os
        os.remove('tmp')

        collection.store.store(
            U, "ler", 'U', collection=collection, eps=eps, rcut=rcut, nmax=nmax, lmax=lmax, metric=metric, n_trees=n_trees, search_k=search_k, **kwargs)

    # Calculate LER
    aid = atoms.get_array('aid')[0]
    result = np.zeros(len(U['clusters']))
    for uid, (center, cluster) in enumerate(U['clusters'].items()):
        result[uid] = np.count_nonzero(np.array([c[0] for c in cluster]) == aid)
    result = result / np.sum(result)
    return result

def test(self, **args):
    return 'test result'
