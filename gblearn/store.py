"""Functions and ResultStore class for storing atomic descriptions
"""
import os
import numpy as np


class ResultStore:
    """

    """

    def __init__(self, location):
        self.root = location
        if not os.path.exists(location):
            os.mkdir(location)

    def check_existing_results(self, descriptor, aid, check_name):
        """ Function to check if correct file structure is in place and if a result file
            exists for these parameters
            Returns:
                bool: True if file results already exits, false if they do not
        """
        path = os.path.join(self.root, descriptor, aid)
        if not os.path.exists(path):
            return False
        for f in os.listdir(path):
            nm, ex = os.path.splitext(f)
            if nm == check_name and os.path.isfile(os.path.join(path, f)):
                return True
        return False

    def store_descriptor(self, result, descriptor, aid, fname, file_ext=None):
        path = os.path.join(self.root, descriptor)
        if not os.path.exists(path):
            os.mkdir(path)
        path = os.path.join(path, aid)
        if not os.path.exists(path):
            os.mkdir(path)
        if isinstance(result, np.ndarray):
            full_path = os.path.join(path, fname)
            np.save(full_path, result)
        else:
            if file_ext is None:
                full_path = full_path + ".dat"
            else:
                full_path = full_path + file_ext
            f = open(full_path, "w+")
            f.write(result)
            f.close()
        return path

    def generate_file_name(self, descriptor, aid, **kwargs):
        """Function to generate file name for storage

        Args:
            descriptor (str): name of descriptor
            aid (str): atoms id 
            **kwargs: arguments for computing descriptor, will be used to generate
                file names

        Returns:
            file name (without extension)

        """
        sep = ""
        name = sep.join([descriptor, "__", aid])
        for key in sorted(kwargs.keys()):
            name = sep.join([name, "___", key, "_", str(kwargs[key])])

        return name

    def get_descriptor(self, aid, descriptor, **args):
        """Function to retrieve descriptor from ResultStore

        Args:
            aid (str):
            descriptor (str):
            **args: refers to whatever arguments that were used to generate the description (which
                correspond to the file name)

        Returns:
            the descriptor, or -1 if no corresponding descriptor is found
        """
        name = self.generate_file_name(descriptor, aid, **args)
        path = os.path.join(self.root, descriptor, aid)
        for f in os.listdir(path):
            nm, ex = os.path.splitext(f)
            if nm == name:
                if ex == ".npy":
                    return np.load(os.path.join(path, f))
                else:
                    return np.loadtxt(os.path.join(path, f))
        return -1
