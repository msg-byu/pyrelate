"""Functions and AtomsCollection class for interacting with collections of Atoms objects
"""
import numpy as np
from tqdm import tqdm
from os import path
import os
from ase import io, Atoms
from gblearn.store import ResultStore


class AtomsCollection(dict):
    """Represents a collection of ASE Atoms objects

    Args :
        name (str) : identifier for this collection, which will be correlated with the ResultStore
            that is created

        fpath (str) : relative file path of where you would like the results to be stored

    Attributes :
        atoms_files (dict) : holds the different ASE atoms objects that make up the collection,
            keys are either user specified (e.g. using regex) or is given by the file name of the
            original data

        r_store (gblearn.store.ResultStore) :

    ''warning'': MAKE SURE TO HAVE UNIQUE COLLECTION NAMES, WILL BE USED FOR LER
    """

    def __init__(self, name):
        super(AtomsCollection, self).__init__()
        self.name = name.lower()

    def _get_aid(self, fpath, comp_rxid, prefix=None):
        """Private function to create the aid for the Atoms object
        """
        extra, fname = path.split(fpath)
        if comp_rxid is not None:
            aid_match = comp_rxid.match(fname)
            if aid_match is None:
                print("Regex found no pattern. Resolving to filename as aid.")
                aid = fname
            else:
                aid = aid_match.group(1)
        else:
            aid = fname

        if prefix is not None:
            aid = prefix.lower() + aid

        return aid

    def read(self, root, Z, f_format=None, rxid=None, prefix=None):
        """Function to read atoms data into ASE Atoms objects and add to Collection

        Args :
            root (str) : relative file path (or list of file paths) to the file, or
                directory of files, where the raw atomic descriptions are located.
            f_format (str) : format of data file. See ase documentation at
                'https://wiki.fysik.dtu.dk/ase/ase/io/io.html'
            Z (int) : atomic number of the elements to be read
            rxid (:obj: str, optional) : regex pattern for extracting the `aid` for each
                Atoms object. Defaults to None. Any files that don't match the regex are
                automatically excluded. The regex should include a named group `(?P<aid>...)`
                so that the id can be extracted correctly. If not specified, the file name is
                used as the `aid`.
            prefix (str): optional prefix for aid

        Example:

           """
        if rxid is not None:
            import re
            comp_rxid = re.compile(rxid)
        try:
            if isinstance(root, list):
                for i in tqdm(range(len(root))):
                    if not isinstance(Z, list):
                        self.read(root[i], Z, f_format, rxid, prefix)
                    else:
                        self.read(root[i], Z[i], f_format, rxid, prefix)
            elif(path.isfile(root)):
                a = io.read(root, format=f_format)
                a.set_atomic_numbers([Z for i in a])
                aid = self._get_aid(root, comp_rxid, prefix)
                self[aid] = a
            elif(path.isdir(root)):
                for afile in os.listdir(root):
                    fpath = os.path.join(root, afile)
                    self.read(fpath, Z, f_format, rxid, prefix)
            else:
                raise ValueError(root)
        except ValueError:
            print("Invalid file path, ", root, " was not read.")

    def describe(self, descriptor, result_store, fcn=None, file_extension=None, **args):
        """Function to call specified description function and store the result

        Args :
            descriptor (str): descriptor to be applied to AtomsCollection, will be used in
                creation of the ResultStore file structure.
            fcn (str): function to apply said description. Built in functions are held in
                descriptors.py, see its documentation for function details.
            atomic_env_specific (bool): Corresponds to if the description yielded by the
                descriptor function is atomic environment specific (True) or collection
                specific (False). Default is true.
            file_extension (str): preferred file extension of stored results. If numpy array,
                extension will be .npy, else the default is .dat
            **args: Parameters associated with the description function specified.

        Returns:

        """
        if fcn is None:
            from gblearn import descriptors
            fcn = getattr(descriptors, descriptor)

        for aid in tqdm(self):
            z = self[aid].get_chemical_symbols()
            fname = result_store.generate_file_name(descriptor, aid, **args)
            exists = result_store.check_existing_results(
                descriptor, aid, fname)
            if not exists:
                result = fcn(self[aid], species=z, **args)
                result_store.store_descriptor(
                    result, descriptor, aid, fname, file_ext=None)
