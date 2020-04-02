"""Functions and Store class for storing atomic descriptions
"""
import os
import numpy as np
import pickle

class Store:
    """Class for efficient storing of description of the AtomsCollection"""

    def __init__(self, location):
        self.root = location
        if not os.path.exists(location):
            os.mkdir(location)

    def _generate_file_name(self, descriptor, idd, **kwargs):
        """Function to generate file name for storage

        Args:
            descriptor (str): name of descriptor
            idd (str): atoms id
            **kwargs (dict): arguments for computing descriptor, will be used to generate
                file names

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

        """ Function to check if correct file structure is in place and if a result file
            exists for these parameters

            Args:
                descriptor (str): name of descriptor
                aid (str): atoms id
                check_name (string): expected filename of results

            Returns:
                bool: True if file results already exits, false if they do not
        """
        fname = self._generate_file_name(descriptor, idd, **kwargs)
        path = os.path.join(self.root, descriptor, idd, fname)
        return os.path.exists(path)

    def store(self, result, descriptor, idd, **kwargs):
        """
        Function to store information into result store

        Args:
            result (pickle): computed result
            descriptor (str): name of descriptor
            idd (str): atoms id
            **kwargs (dict): arguments for computing descriptor, will be
                used to generate file names

        """
        fname = self._generate_file_name(descriptor, idd, **kwargs)
        path = os.path.join(self.root, descriptor, idd)
        os.makedirs(path, exist_ok=True)

        path = os.path.join(path, fname)
        with open(path, 'wb') as f:
            pickle.dump(result, f)

    def get(self, descriptor, idd, **kwargs):
        """Function to retrieve information from the Store

        Args
            descriptor (str): name of descriptor
            idd (str): atoms id
            **kwargs (dict): refers to whatever arguments that were used to generate
                the description (which correspond to the file name)
        """
        fname = self._generate_file_name(descriptor, idd, **kwargs)
        path = os.path.join(self.root, descriptor, idd, fname)

        result = None
        try:
            with open(path, 'rb') as f:
                result = pickle.load(f)
        except FileNotFoundError:
            pass

        return result
