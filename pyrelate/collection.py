"Functions and AtomsCollection class for interacting with collections of ASE Atoms objects"""
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

    def __init__(self, name, store=None, data=None):
        """Initializer which calls dict's and Store's initializers."""
        super(AtomsCollection, self).__init__()
        self.name = name.lower()
        #TODO make trim and pad more general (have FIRM for entire collection, apply to new atoms read in, cannot trim
        # to be wider than it is, when you do "subset" the trim/pad values must be applied, update trim function)
        # self.trim = trim
        # self.pad = pad

        if type(store) == Store:
            self.store = store
        elif type(store) == str:
            self.store = Store(store)

        if data is not None:
            self.update(data) #allows you to initialize a collection beginning with another collection (dictionary)

    def __str__(self):
        """String representation of the AtomsCollection object (name of collection)."""
        return self.name

    def subset(self, aids, name=None, store=None):
        """Return an AtomsCollection containing the specified subset of the original collection.

        Parameters:
            aids (list): List containing the strings of all the aids to include in the new collection.
            name (string): Name of the new collection, default is the name of the original collection.
            store (string): String representing the location of the store, default (None) sets the current
                store as the store for the new collection, which is probably what you want most of the time.

        Returns:
            pyrelate.AtomsCollection
        """
        if name is None:
            name = self.name
        if store is None:
            store = self.store
        data = {aid: self[aid] for aid in aids}
        return AtomsCollection(name, store=store, data=data)

    def _read_aid(self, fpath, comp_rxid, prefix=None):
        """Private function to read the aid for the Atoms object from filename.

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

    def read(self, root, Z, f_format=None, rxid=None, prefix=None):
        """Function to read atoms data into ASE Atoms objects and add to AtomsCollection.

        Utilizes functionality in ASE to read in atomic data.

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

                #NEW - initialize mask to all ones (keep all)
                a.new_array("mask", np.array([1 for i in range(len(a))]), dtype="int")
            elif(path.isdir(root)):
                for afile in os.listdir(root):
                    fpath = os.path.join(root, afile)
                    self.read(fpath, Z, f_format, rxid, prefix)
            else:
                raise ValueError(root)
        except ValueError:
            print("Invalid file path,", root, "was not read.")

    def trim(self, trim, dim, pad=True):
        """Trims off excess atoms and indicates padding (specified in a mask).
        #FIXME store trim and pad values for the entire collection

        Padding may want to be included so that atoms have a full atomic environments at the edge, some descriptors that perform better with full atomic environments. The "mask" array is attached to the Atoms object, with 1 indicating atoms to be kept in the final description, and 0 to indicate padding atoms. You may get the mask by calling 'my_col.get_array("mask")'.

        Parameters:
            trim (float or int): value (in Angstroms) of the furthest atoms that you want included in the calculated results
            dim (int): what dimension the grain boundary is in, 0 for x, 1 for y, 2 for z
            pad (boolean or float or int): 'True' (default) gives a padding value equal to that of the trim, 'False' gives no padding, or specify the amount of padding wanted (in Angstroms).

        .. WARNING:: The user must keep track of what trim and pad parameters were used, as the results stored will not indicate if they have been trimmed or not

        Example:
            .. code-block:: python

                my_col.trim(trim=4, dim=0)
                my_col.trim(4, 0, pad=False)
        """

        if pad is True:
            pad = trim
        elif pad is False:
            pad = 0

        if not (isinstance(trim, int) or isinstance(trim, float)):
            raise TypeError("Trim should be int or float type")
        if not (isinstance(pad, int) or isinstance(pad, float)):
            raise TypeError("Pad should be int, float, or boolean type")
        if not isinstance(dim, list):
            dim = [dim for i in range(len(self))]

        slice_width = trim + pad
        for idx, aid in enumerate(tqdm(self)):
            atoms = self[aid]
            d = dim[idx]
            if d not in [0,1,2]:
                raise TypeError("Dimension should equal 0, 1, or 2")

            pos = atoms.get_positions()[:,d]
            #TODO verify gbcenter = 0
            gbcenter = 0
            #delete atoms outside of trim and pad
            del atoms[[atom.index for atom in atoms if atom.position[d] < (gbcenter - slice_width) or atom.position[d] > (gbcenter + slice_width) ]]

            #update mask -- 1 for inside trim, and 0 for pad
            mask = np.array([atom.position[d] > (gbcenter - trim) for atom in atoms])*np.array([atom.position[d] < (gbcenter + trim) for atom in atoms])*1
            atoms.set_array("mask", mask)

    def describe(self, descriptor, aid=None, fcn=None, override=False, **desc_args):
        """Function to calculate and store atomic description.

        User can specify a descriptor function to be used, or use those in descriptors.py. When there is a padding associated with the Atoms object, the padding atoms are deleted from the final description before being stored.

        Parameters:
            descriptor (str): descriptor to be applied.
            aid (iterable or string): Atoms ID's (aid) of atomic systems to be described. Can pass in single aid string, or an iterable
            of aids. Defaults to None. When None, all ASE Atoms objects in AtomsCollection are described.
            fcn: function to apply said description. Defaults to none. When none, built in functions in descriptors.py are used.
            override (bool): if True, descriptor will override any matching results in the store. Defaults to False.
            desc_args (dict): Parameters associated with the description function specified. See documentation in descriptors.py for function details and parameters.

        Examples:
            #FIXME old examples
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
                    "res_needed": "my_soap",
                    "soap_fcn": soap,
                    "eps": 0.1,
                    "dissim_args":{"gamma":4000},
                    **soap_args
                }
                my_col.describe("ler", **ler_args)

        """

        if fcn is None:
            from pyrelate import descriptors
            fcn = getattr(descriptors, descriptor)

        if aid is None:
            to_calculate = self.aids()
        else:
            to_calculate = [aid] if type(aid) is str else aid

        for aid in tqdm(to_calculate):
            if aid not in self.aids():
                raise ValueError(f"{aid} is not a valid atoms ID.")

            exists = self.store.check_exists("Descriptions", descriptor, aid, **desc_args)
            if not exists or override:
                result = fcn(self[aid], **desc_args)

            if result is not None:
                if len(result) > np.count_nonzero(self[aid].get_array( "mask")):
                    to_delete = np.logical_not(self[aid].get_array("mask"))
                    result = np.delete(result, to_delete, axis=0)

                #FIXME store trim/pad data in info dict
                info = {} #"trim":None, "pad":None}

                self.store.store_description(
                    result, info, aid, descriptor,  **desc_args)
        # self.clear("temp")

    def process(self, method, based_on, fcn=None, override=None, **kwargs):
        """Calculate and store collection specific results

        Parameters:
            method (str): string indicating the name of the descriptor to be applied
            based_on (str, dict): tuple holding the information necessary to get previously calculated results for use in processing.
            fcn (str): function to apply said processing method. Defaults to none. When none, built in functions in descriptors.py are used.
            override (bool): if True, results will override any matching results in the store. Defaults to False.
            kwargs (dict): Parameters associated with the processing function specified. See documentation in descriptors.py for function details and parameters.

        Examples:

        """

        if fcn is None:
            from pyrelate import descriptors
            fcn = getattr(descriptors, method)

        # TODO: if soap not calculated, calculate
        # would require the function used to be passed in somehow...
        # desc = based_on[0]
        # desc_args = based_on[1]
        # for aid in self.aids():
        #     if not self.store.check_exists("Descriptions", aid, desc, **desc_args):
        #         self.describe(desc, aid=aid, **desc_args)

        exists = self.store.check_exists("Collections", self.name, method, based_on, **kwargs)
        if not exists or override:
            print("Calculate and store")
            result, info = fcn(self, method, based_on, **kwargs)
            self.store.store_collection_result(result, info, method, self.name, based_on, **kwargs)

        return self.store.get_collection_result(method, self.name, based_on, **kwargs) #return result and info dict


    def clear(self, descriptor=None, idd=None, **kwargs):
        '''Function to delete specified results from Store.

        Functionality includes:

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
        #FIXME rewrite how change works, and it is implemented differently in store.py
        # has_kwargs = len(kwargs) != 0
        # if descriptor is not None:
        #     if has_kwargs:
        #         if idd is not None:
        #             self.store.clear(descriptor, idd, **kwargs)
        #         else:
        #             idds = self.aids()
        #             self.store.clear(descriptor, idds, **kwargs)
        #     else:
        #         self.store.clear_descriptor(descriptor)
        # else:
        #     self.store.clear_all()
        pass

    def get_description(self, idd, descriptor, **desc_args):
        """Wrapper function to retrieve descriptor results from the store"""
        return self.store.get_description(idd, descriptor, **desc_args)

    def get_collection_result(self, method, collection_name, based_on, **method_args):
        """Wrapper function to retrieve collection specific results from the store"""
        return self.store.get_collection_result(method, collection_name, based_on, **method_args)

    def aids(self):
        '''Returns sorted list of atom ID's (aids) in collection'''
        a = list(self)
        a.sort()
        return a
