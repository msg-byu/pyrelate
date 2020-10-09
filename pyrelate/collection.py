"""Functions and AtomsCollection class for interacting with collections of ASE Atoms objects"""
import numpy as np
from tqdm import tqdm
from os import path
import os
from ase import io, Atoms
from pyrelate.store import Store


class AtomsCollection(dict):
    """Represents a collection of ASE Atoms objects

    :Attributes:
        - self (dict): inherits from dictionary
        - name (str) : identifier for this collection
        - store (Store) : store to hold all the results and other information

    .. WARNING:: Make sure to have unique collection names, will be used for LER

    """

    def __init__(self, name, store_path=None):
        """Initializer which calls dict's and Store's initializers"""
        super(AtomsCollection, self).__init__()
        self.name = name.lower()
        self.store = Store(store_path)

    def __str__(self):
        """String representation of the AtomsCollection object (name of collection)"""
        return self.name

    def _read_aid(self, fpath, comp_rxid, prefix=None):
        """Private function to read the aid for the Atoms object from filename

        Parameters:
            fpath (str): file path to the atomic information to be read in
            comp_rxid(_sre.SRE_Pattern): pre-compiled regex parser to extract desired aid from file name. If none found, default aid will be the file name.
            prefix (str): otional prefix for aid to be generated (will be made lowercase)

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
            aid = prefix.lower() + "_" + aid

        return aid

    def _descriptor_needs_store(self, fcn):
        from inspect import signature
        try:
            signature(fcn).parameters['store']
            return True
        except:
            return False

    def read(self, root, Z, f_format=None, rxid=None, prefix=None):
        """Function to read atoms data into ASE Atoms objects and add to AtomsCollection

        Parameters:
            root (str) : relative file path (or list of file paths) to the file, or directory of files, where the raw atomic descriptions are located.
            Z (int) : atomic number of the elements to be read
            f_format (str) : format of data file. Defaults to None. See ASE's documentation at 'https://wiki.fysik.dtu.dk/ase/ase/io/io.html'
            rxid (:obj: str, optional) : regex pattern for extracting the `aid` for each Atoms object. Defaults to None. The regex should include a named group `(?P<aid>...)` so that the id can be extracted correctly.  If any files don't match the regex or if it is not specified, the file name is used as the `aid`.
            prefix (str): optional prefix for aid. Defaults to none.

        Example:
            .. code-block:: python

                my_col.read(["/Ni/ni.p454.out", "/Ni/ni.p453.out"], 28, "lammps-dump-text", rxid=r'ni.p(?P<aid>\d+).out', prefix="Nickel")
                my_col.read("/Ni/", 28, "lammps-dump-text", rxid=r'ni.p(?P<aid>\d+).out', prefix="Nickel")

        """
        # TODO add functionality to pass a collection into read
        comp_rxid = None
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
                # TODO generalize for reading multi elemental data
                a0 = io.read(root, format=f_format)
                a = a0.copy()
                # delete end blocks
                del a[[atom.index for atom in a if atom.number ==
                       4 or atom.number == 5]]
                a.new_array('type', a.get_array(
                    'numbers', copy=True), dtype=int)
                a.set_atomic_numbers([Z for i in a])
                # TODO aid is stored as an array in the Atoms object, ideally want a single property for the Atoms object
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

    def describe(self, descriptor, fcn=None, override=False, **kwargs):
        """Function to call specified description function and store the result

        Parameters:
            descriptor (str): descriptor to be applied to AtomsCollection.
            fcn (str): function to apply said description. Defaults to none. When none, built in functions in descriptors.py are used.
            override (bool): if True, descriptor will override any matching results in the store. Defaults to False.
            kwargs (dict): Parameters associated with the description function specified. See documentation in descriptors.py for function details and parameters.

        Example:
            .. code-block:: python
            
                soap_args = {
                    "rcut" : 5,
                    "lmax" : 9,
                    "nmax" : 9
                }
                my_col.describe("soap", **soap_args)
                my_col.describe("my_soap", fcn=soap, **soap_args)
                my_col.describe("asr", res_needed="my_soap", **soap_args)

                ler_args = {
                    "collection": my_col,
                    "eps": 0.025,
                    "res_needed": "soap",
                    **soap_args
                }
                my_col.describe("ler", **ler_args)

        """

        if fcn is None:
            from pyrelate import descriptors
            fcn = getattr(descriptors, descriptor)

        for aid in tqdm(self):
            exists = self.store.check_exists(
                descriptor, aid, **kwargs)
            if not exists or override:
                if self._descriptor_needs_store(fcn):
                    result = fcn(self[aid], self.store, **kwargs)
                else:
                    result = fcn(self[aid], **kwargs)

                if result is not None:
                    self.store.store(
                        result, descriptor, aid, **kwargs)
        self.clear("temp")

    def clear(self, descriptor=None, idd=None, **kwargs):
        '''Function to delete specified results from Store. You can:

            - remove a result for a specific Atoms object

            - remove specific results for all Atoms objects in the collection

            - remove all results for a certain type of descriptor, and

            - remove all results in the store

        Parameters:
            descriptor (str): descriptor to be applied to AtomsCollection.
            idd (str): Idd of Atoms object who's results you want. Defaults to None, in which case corresponding results for all ASE Atoms objects will be returned as a dictionary, with the aid as key.
            kwargs (dict): Parameters associated with the description function specified.

        Example:
            .. code-block:: python

                my_col.clear("soap", "aid1", rcut=5.0, nmax=9, lmax=9) #clears single SOAP result with given parameters
                my_col.clear("soap", rcut=5.0, nmax=9, lmax=9) #clears SOAP results for whole collection with given parameters
                my_col.clear("soap") #clears all soap results from store
                my_col.clear() #clears all results from store

        '''
        has_kwargs = len(kwargs) != 0
        if descriptor is not None:
            if has_kwargs:
                if idd is not None:
                    self.store.clear(descriptor, idd, **kwargs)
                else:
                    idds = self.aids()
                    self.store.clear(descriptor, idds, **kwargs)
            else:
                self.store.clear_descriptor(descriptor)
        else:
            self.store.clear_all()

    def get(self, descriptor, idd=None, **kwargs):
        '''Shell function to call Store's get method

        Parameters:
            descriptor (str): descriptor to be applied to AtomsCollection.
            idd (str): Idd of Atoms object who's results you want. Defaults to None, in which case corresponding results for all ASE Atoms objects will be returned as a dictionary, with the aid as key.
            kwargs (dict): Parameters associated with the description function specified.

        Example:
            .. code-block:: python

                my_col.get("soap", 'aid1', rcut=5.0, nmax=9, lmax=9)
                my_col.get("asr", res_needed="my_soap", rcut=5.0, nmax=9, lmax=9)

        '''
        if idd is None:
            idd = self.aids()

        return self.store.get(descriptor, idd, **kwargs)

    def aids(self):
        '''Returns sorted list of atom ID's (aids) in collection'''
        a = list(self)
        a.sort()
        return a
