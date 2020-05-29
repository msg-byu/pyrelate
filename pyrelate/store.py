"""Functions and Store class for storing atomic descriptions
"""
import os
import numpy as np
import pickle


class Store:
    """Class for efficient storing of description of the AtomsCollection"""

    def __init__(self, location=None):
        """Default initiation creates a Store entitled 'store' in current directory"""
        if location is None:
            self.root = os.path.join(os.getcwd(), "store")
        else:
            self.root = location
        if not os.path.exists(self.root):
            os.mkdir(self.root)

    def __str__(self):
        """Returns relative path of the store location for the string representation of Store object"""
        return os.path.relpath(self.root)

    def _generate_file_name(self, descriptor, idd, **kwargs):
        """Function to generate file name for storage

        Parameters:
            descriptor (str): name of descriptor
            idd (str): atoms id
            kwargs (dict): arguments for computing descriptor, will be used to generate file names

        Returns:
            string: file name with extension '.pkl'

        """
        sep = ""
        name = sep.join([descriptor, "__", idd])
        for key in sorted(kwargs.keys()):
            name = sep.join([name, "___", key, "_", str(kwargs[key])])
        name += '.pkl'

        return name

    def check_exists(self, descriptor, idd, **kwargs):
        """ Function to check if correct file structure is in place and if a result file exists for these parameters

        Parameters:
            descriptor (str): name of descriptor.
            idd (str): id of item to check for.
            kwargs (dict): Parameters associated with the description function.

        Returns:
            bool: True if file results already exists, false if they do not
        """
        fname = self._generate_file_name(descriptor, idd, **kwargs)
        path = os.path.join(self.root, descriptor, idd, fname)
        return os.path.exists(path)

    def store(self, result, descriptor, idd, **kwargs):
        """
        Function to store information into result store

        Parameters:
            result (pickle): computed result
            descriptor (str): name of descriptor
            idd (str): atoms id
            kwargs (dict): arguments for computing descriptor, will be used to generate file names

        """
        fname = self._generate_file_name(descriptor, idd, **kwargs)
        path = os.path.join(self.root, descriptor, idd)
        os.makedirs(path, exist_ok=True)

        path = os.path.join(path, fname)
        with open(path, 'wb') as f:
            pickle.dump(result, f)

    def _get_file(self, descriptor, idd, **kwargs):
        """Function to retrieve information from a file in the Store

        Parameters:
            descriptor (str): name of descriptor
            idd (str): atoms id
            kwargs (dict): refers to whatever arguments that were used to generate the description (which correspond to the file name)
        """

        fname = self._generate_file_name(descriptor, idd, **kwargs)
        path = os.path.join(self.root, descriptor, idd, fname)

        result = None
        try:
            with open(path, 'rb') as f:
                result = pickle.load(f)
        except FileNotFoundError:
            pass
        except Exception as exception:
            print("%s when loading file %s, consider deleting result and recomputing" % (type(exception).__name__, fname))

        return result

    def get(self, descriptor, idd, **kwargs):
        """Function to retrieve information from the Store

        Parameters:
            descriptor (str): name of descriptor
            idd (str **or** list(str)): atoms id (or list of atom id's)
            kwargs (dict): refers to whatever arguments that were used to generate the description (which correspond to the file name)
        """
        if type(idd) is list:
            result ={}
            for i in idd:
                result[i] = self._get_file(descriptor, i, **kwargs)
        else:
            result = self._get_file(descriptor, idd, **kwargs)

        return result
