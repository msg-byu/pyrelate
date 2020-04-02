"""Functions and AtomsCollection class for interacting with collections of Atoms objects
"""
import numpy as np
from tqdm import tqdm
from os import path
import os
from ase import io, Atoms
from gblearn.store import Store


class AtomsCollection(dict):
    """Represents a collection of ASE Atoms objects

    Attributes :
        name (str) : identifier for this collection
        self (dict): inherits from dictionary

    ..warning:: MAKE SURE TO HAVE UNIQUE COLLECTION NAMES, WILL BE USED FOR LER
    """

    def __init__(self, name, store_path):
        """initializer which calls dicts initializer"""
        super(AtomsCollection, self).__init__()
        self.name = name.lower()
        self.store = Store(store_path)

    def __str__(self):
        return self.name

    def _read_aid(self, fpath, comp_rxid, prefix=None):
        """Private function to create the aid for the Atoms object

        Args:
            fpath (str): file path to the file holding the atoms information to
                be read in. File name will be used in aid creation.
            comp_rxid(_sre.SRE_Pattern): pre-compiled regex parser to extract desired
                aid from file name. If none found, default aid will be the file name.
            prefix (str): otional prefix for aid to be generated

        Returns:
            aid (str): atoms id, will be used as key for the corresponding Atoms object
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
            Z (int) : atomic number of the elements to be read
            f_format (str) : format of data file. Defaults to None. See ase documentation at
                'https://wiki.fysik.dtu.dk/ase/ase/io/io.html'
            rxid (:obj: str, optional) : regex pattern for extracting the `aid` for each
                Atoms object. Defaults to None. Any files that don't match the regex are
                automatically excluded. The regex should include a named group `(?P<aid>...)`
                so that the id can be extracted correctly. If not specified, the file name is
                used as the `aid`.
            prefix (str): optional prefix for aid

        Example:
            c.read(["../homer/ni.p454.out", "../homer/ni.p453.out"], 28,
                "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out',
                prefix="Homer")

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
                #FIXME store aid as array of atoms object, do we want that?
                aid = self._read_aid(root, comp_rxid, prefix)
                a.new_array("aid", [aid for i in a], dtype="str")
                self[aid] = a
            elif(path.isdir(root)):
                for afile in os.listdir(root):
                    fpath = os.path.join(root, afile)
                    self.read(fpath, Z, f_format, rxid, prefix)
            else:
                raise ValueError(root)
        except ValueError:
            print("Invalid file path,", root, "was not read.")

    def describe(self, descriptor, fcn=None, needs_store=False, **kwargs):
        #FIXME check if result store is not none
        """Function to call specified description function and store the result

        Args :
            descriptor (str): descriptor to be applied to AtomsCollection, will be used in
                creation of the ResultStore file structure.
            store (gblearn.store.ResultStore): Object to facilitate storage of
                the results of describe
            fcn (str): function to apply said description. Built in functions are held in
                descriptors.py, see its documentation for function details.
            **kwargs (dict): Parameters associated with the description function specified.

        Returns:
            None: everything will be stored in the specified location in Result Store
        Example:
            #FIXME: change example to updated version
            rs = ResultStore("../store")
            c.describe("soap", rs,  rcut=5.0, nmax=9, lmax=9)
        """

        if fcn is None:
            from gblearn import descriptors
            fcn = getattr(descriptors, descriptor)

        for aid in tqdm(self):
            exists = self.store.check_exists(
                descriptor, aid, **kwargs)
            if not exists:
                if needs_store:
                    result = fcn(self[aid], self.store, **kwargs)
                else:
                    result = fcn(self[aid], **kwargs)
                self.store.store(
                    result, descriptor, aid, **kwargs)
