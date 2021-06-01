import numpy as np
'''Built-in descriptors for use with AtomsCollection's describe function.

Some guidelines for writing your own descriptor function
    - Some descriptions (including built in ASR and LER) require previously computed results, if you want to write your own descriptor function that accesses the store, one of your included parameters must be 'store'. 'Res_needed' is also a suggested parameter to include, it indicates which results to use when computing the new descriptor.
'''


def soap(atoms, rcut, nmax, lmax, **kwargs):
    """Smooth Overlap of Atomic Positions-- pycsoap implementation

    Parameters:
        atoms ('ase.atoms.Atoms'): ASE atoms object to perform description on
        nmax (int): bandwidth limits for the SOAP descriptor radial basis functions.
        lmax (int): bandwidth limits for the SOAP descriptor spherical harmonics.
        rcut (float): local environment finite cutoff parameter.
        kwargs (dict): Parameters associated with the description function

    To see `pycsoap's documentation <https://pypi.org/project/pycsoap/>`_.
    """
    from pycsoap.soaplite import SOAP
    soap_desc = SOAP(atomic_numbers=atoms.numbers, rcut=rcut,
                     nmax=nmax, lmax=lmax, **kwargs)
    P = soap_desc.create(atoms)
    return P


def asr(collection, based_on, norm_asr=False):
    """Average SOAP representation: average vectors from SOAP matrix into a single vector

    Parameters:
        norm_asr (bool): Normalize ASR vector. Default is False, not normalized.
    """
    # for each gb, average and add to matrix, and return
    asr_matrix = []
    magnitude = 1
    for aid in collection.aids():
        res = collection.get_description(aid, based_on[0], **based_on[1])
        asr_row = np.average(res, axis=0)
        if norm_asr is True:
            magnitude = np.linalg.norm(asr_row)
        asr_matrix.append(asr_row / magnitude)
    return np.array(asr_matrix)


def sum(collection, based_on):
    """Sum all rows of a descriptor matrix into a single vector

    Parameters:

    """
    sum_matrix = []
    for aid in collection.aids():
        res = collection.get_description(aid, based_on[0], **based_on[1])
        sum_row = np.sum(res, axis=0)
        sum_matrix.append(sum_row)
    return np.array(sum_matrix)


def gaussian_dissimilarity(lae1, lae2, gamma):
    """Gaussian (RBF) dissimilarity metric (used with LER) that varies between 0 and 1.

    This metric is given by :math:`1 - \exp{(-\gamma*||x - x'||^2)}`. A dissimilarity value closer to 0 means the LAEs are more similar, while closer to 1 means more dissimilar.

    Parameters:
        lae1 (np.array): first lae to compare
        lae2 (np.array): second lae to compare
        gamma (float): gamma value
    """
    diff = lae1 - lae2
    gker = 1 - np.exp((-gamma * (np.linalg.norm(diff)**2)))
    return gker

# def ler(atoms, store, collection, eps, res_needed='soap', dissimilarity=gaussian_dissimilarity, dissim_args={}, soap_fcn=None, seed=None, metric='euclidean', n_trees=10, search_k=-1, **soapargs):
#     '''Local Environment Representation
#
#         Parameters:
#             atoms ('ase.atoms.Atoms'): ASE atoms object to perform description on
#             store (pyrelate.store.Store) : store to access previously computed SOAP matrix (automatically passed in with AtomsCollection.describe())
#             collection (pyrelate.collection.AtomsCollection): LER is collection specific, needed for computation
#             eps (float): epsilon value indicating that any LAE's outside this value are considered dissimilar. Descriptor and dissimilarity metric specific.
#             res_needed (str): name of description whose results are needed to compute the LER. Defaults to 'soap'.
#             dissimilarity (function): dissimilarity metric function that has the first 2 parameters as the 2 LAE's to be compared.
#             dissim_args (dict): dictionary with any additional hyperparameter arguments for the given dissimilarity metric.
#             soap_fcn (function): optional parameter for a function to compute SOAP matrix for the element's perfect crystal on the fly. Defaults to None. When None, 'soap' function in descriptors.py will be used.
#             seed(np.ndarray or list): perfect seed for the element being considered. Defaults to None. When None, seed will be generated on the fly.
#             metric(str): For approximate nearest neighbor calculation. See Annoy's documentation for more details.
#             n_trees(int): For approximate nearest neighbor calculation. See Annoy's documentation for more details.
#             search_k(int): For approximate nearest neighbor calculation. See Annoy's documentation for more details.
#             soapargs (dict): Parameters associated with the SOAP description being used.
#         `Annoy's documentation <https://github.com/spotify/annoy>`_.
#     '''
#
#     U = store.get(
#         "temp", 'U_ler', collection=collection, eps=eps, dissim_args=dissim_args, res_needed=res_needed, soap_fcn=soap_fcn, metric=metric, n_trees=n_trees, search_k=search_k, **soapargs)  # add seed? or hash for seed?
#
#     if U is None:
#         from collections import OrderedDict
#         U = {
#             'centers': OrderedDict(),
#             'clusters': {}
#         }
#
#         # Part 1: Clustering
#         import pyrelate.elements as elements
#         if seed is None:
#             seed = elements.seed(list(collection.values())[0].get_chemical_symbols()[0], soap_fcn, **soapargs)
#         U['centers'][('0', 0)] = seed
#         #print("Part 1: Clustering (for each LAE in each GB)")
#         for aid in collection:
#             for lae_num, lae in enumerate(store.get(res_needed, aid, **soapargs)):
#                 if lae is None:
#                     raise RuntimeError(
#                         "LER requires SOAP to be generated first")
#                 for unique in U['centers'].values():
#                     dist = dissimilarity(unique, lae, **dissim_args)
#                     if dist < eps:
#                         break
#                 else:
#                     U['centers'][(aid, lae_num)] = np.copy(lae)
#         # TODO sort the unique bins
#
#         # Part 2: Classifying
#         from annoy import AnnoyIndex
#         dim = len(seed)
#         a = AnnoyIndex(dim, metric)
#         for i, unique in enumerate(U['centers'].values()):
#             a.add_item(i, unique)
#         a.build(n_trees)
#         a.save('tmp')  # TODO find better way to save temporary file
#         for aid in collection:
#             for lae_num, lae in enumerate(store.get(res_needed, aid, **soapargs)):
#                 nni = a.get_nns_by_vector(
#                     lae, 1, search_k=search_k, include_distances=False)[0]
#                 cluster = list(U['centers'].keys())[nni]
#                 # TODO when we sort just fill the empty lists first
#                 if cluster not in U['clusters']:
#                     U['clusters'][cluster] = []
#                 U['clusters'][cluster].append((aid, lae_num))
#         import os
#         os.remove('tmp')
#
#         store.store(
#             U, "temp", 'U_ler', collection=collection, eps=eps, dissim_args=dissim_args, res_needed=res_needed, soap_fcn=soap_fcn, metric=metric, n_trees=n_trees, search_k=search_k, **soapargs)
#
#     # Calculate LER
#     aid = atoms.get_array('aid')[0]
#     result = np.zeros(len(U['clusters']))
#     for uid, (center, cluster) in enumerate(U['clusters'].items()):
#         result[uid] = np.count_nonzero(
#             np.array([c[0] for c in cluster]) == aid)
#     result = result / np.sum(result)
#     return result


def ler(collection, based_on, eps, dissimilarity=gaussian_dissimilarity, dissim_args={}, soap_fcn=None, seed=None, metric='euclidean', n_trees=10, search_k=-1, **kwargs):
    '''Local Environment Representation

        Parameters:
            atoms ('ase.atoms.Atoms'): ASE atoms object to perform description on
            store (pyrelate.store.Store) : store to access previously computed SOAP matrix (automatically passed in with AtomsCollection.describe())
            collection (pyrelate.collection.AtomsCollection): LER is collection specific, needed for computation
            eps (float): epsilon value indicating that any LAE's outside this value are considered dissimilar. Descriptor and dissimilarity metric specific.
            res_needed (str): name of description whose results are needed to compute the LER. Defaults to 'soap'.
            dissimilarity (function): dissimilarity metric function that has the first 2 parameters as the 2 LAE's to be compared.
            dissim_args (dict): dictionary with any additional hyperparameter arguments for the given dissimilarity metric.
            soap_fcn (function): optional parameter for a function to compute SOAP matrix for the element's perfect crystal on the fly. Defaults to None. When None, 'soap' function in descriptors.py will be used.
            seed(np.ndarray or list): perfect seed for the element being considered. Defaults to None. When None, seed will be generated on the fly.
            metric(str): For approximate nearest neighbor calculation. See Annoy's documentation for more details.
            n_trees(int): For approximate nearest neighbor calculation. See Annoy's documentation for more details.
            search_k(int): For approximate nearest neighbor calculation. See Annoy's documentation for more details.
            soapargs (dict): Parameters associated with the SOAP description being used.
            `Annoy's documentation <https://github.com/spotify/annoy>`_.
    '''

    from collections import OrderedDict
    U = {
        'centers': OrderedDict(),
        'clusters': {}
    }

    # Part 1: Clustering
    import pyrelate.elements as elements
    if seed is None:
        seed = elements.seed(list(collection.values())[0].get_chemical_symbols()[0], soap_fcn, **based_on[1])
    U['centers'][('0', 0)] = seed
    # print("Part 1: Clustering (for each LAE in each GB)")
    for aid in collection:
        soap = collection.get_description(aid, based_on[0], **based_on[1])
        for lae_num, lae in enumerate(soap):
            if lae is None:
                raise RuntimeError(
                    "LER requires SOAP to be generated first")
            for unique in U['centers'].values():
                dist = dissimilarity(unique, lae, **dissim_args)
                if dist < eps:
                    break
            else:
                U['centers'][(aid, lae_num)] = np.copy(lae)
    # TODO sort the unique bins

    # Part 2: Classifying
    from annoy import AnnoyIndex
    dim = len(seed)
    a = AnnoyIndex(dim, metric)
    for i, unique in enumerate(U['centers'].values()):
        a.add_item(i, unique)
    a.build(n_trees)
    a.save('tmp')  # TODO find better way to save temporary file
    for aid in collection:
        soap = collection.get_description(aid, based_on[0], **based_on[1])
        for lae_num, lae in enumerate(soap):
            nni = a.get_nns_by_vector(
                lae, 1, search_k=search_k, include_distances=False)[0]
            cluster = list(U['centers'].keys())[nni]
            # TODO when we sort just fill the empty lists first
            if cluster not in U['clusters']:
                U['clusters'][cluster] = []
            U['clusters'][cluster].append((aid, lae_num))
    import os
    os.remove('tmp')

    # Calculate LER
    ler_matrix = []
    for aid in collection.aids():
        result = np.zeros(len(U['clusters']))
        for uid, (center, cluster) in enumerate(U['clusters'].items()):
            result[uid] = np.count_nonzero(
                np.array([c[0] for c in cluster]) == aid)
        result = result / np.sum(result)
        ler_matrix.append(result)

    info = {
        "num_clusters": len(U['centers']),
        "cluster_centers": U['centers']
    }
    return np.array(ler_matrix), info
