"""Functions and Collection class for interacting with collections of Atoms objects
"""
#import numpy as np
from tqdm import tqdm
from os import path
import os
from ase import io, Atoms
from gblearn.store import ResultStore

class Collection:
    #inherit from a dictionary, (look in gblearn to double check), can get rid of atoms_files
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

    def __init__(self, name, fpath):
        super(Collection, self).__init__()
        self.atoms_files = {}
        self.name = name.lower()
        self._result_store = ResultStore(fpath)

    def _get_idd(self, fpath,comp_rxid, prefix=None):
        """Private function to create the idd for the Atoms object
        """
        extra,fname = path.split(fpath)
        if comp_rxid is not None:
            idd_match = comp_rxid.match(fname)
            if idd_match is None:
                print("Regex found no pattern. Resolving to filename as idd.")
                idd = fname
            else:
                idd = idd_match.group(1)
        else: 
            idd=fname
        
        if prefix is not None:
            idd= prefix.lower() + idd

        return idd

    def read(self, root, Z, f_format=None, rxid=None, prefix=None):
        """Function to read atoms data into ASE Atoms objects and add to Collection

        Args :
            root (str) : relative file path (or list of file paths) to the file, or 
                directory of files, where the raw atomic descriptions are located.
            f_format (str) : format of data file. See ase documentation at
                'https://wiki.fysik.dtu.dk/ase/ase/io/io.html'
            Z (int) : atomic number of the elements to be read
            rxid (:obj: str, optional) : regex pattern for extracting the `idd` for each 
                Atoms object. Defaults to None. Any files that don't match the regex are 
                automatically excluded. The regex should include a named group `(?P<idd>...)` 
                so that the id can be extracted correctly. If not specified, the file name is 
                used as the `idd`.
            prefix (str): optional prefix for idd

        Example:

           """
        if rxid is not None:
            import re
            comp_rxid = re.compile(rxid)
        try:
            if isinstance(root, list):
                for i in range(len(root)):
                    if not isinstance(Z, list):
                        self.read(root[i], f_format, Z, rxid, prefix)
                    else:
                        self.read(root[i], f_format,Z[i], rxid, prefix) 
            elif(path.isfile(root)):
                a=io.read(root, format=f_format)
                a.set_atomic_numbers([Z for i in a])
                idd = self._get_idd(root, comp_rxid, prefix)
                self.atoms_files[idd] = a
            elif(path.isdir(root)):
                for afile in os.listdir(root):
                    fpath = os.path.join(root, afile)
                    self.read(fpath, f_format, Z, rxid, prefix)
            else:
                raise ValueError(root)
        except ValueError:
            print("Invalid file path, " , root, " was not read.")

    def describe(self, descriptor, fcn=None, atomic_env_specific=True, file_extension=None , **args):
        """Function to call specified description function and store the result
        
        Args :
            descriptor (str): descriptor to be applied to Collection, will be used in 
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
            #TODO: work on it
        
        for idd in self.atoms_files:
            z = self.atoms_files[idd].get_chemical_symbols() #.get_atomic_numbers()
            fname = self._result_store.generate_file_name(descriptor, idd,z[0], **args)
            exists, fpath =  self._result_store.check_existing_results(descriptor, idd, fname, atomic_env_specific)
            if not exists:
                result = fcn(self.atoms_files[idd], species=z, **args)
                self._result_store.store_to_file(result, fpath, fname, file_ext=None)

#Questions
#what else do I need before I push to gblearn? (unit tests, travis CI)
#can you look over everything pls?

#TODO
#better error catching
#have user input where to get module that corresponds to their own describe function

#Documentation--docstring, look at correct syntax, compare to other documentation, 
#link to other documentation

#Get descriptors working more generally
#add tqdm bar for loading GB's??
#go through getting started and follow structure of developing scientific code

