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
    
    def check_existing_results(self, descriptor, idd, check_name, atomic_env_specific):
        """ Function to check if correct file structure is in place and if a result file
            exists for these parameters
            Returns:
                bool: True if file results already exits, false if they do not
        """    
        path = os.path.join(self.root, descriptor)
        if not os.path.exists(path): #if descriptor dir exists
            os.mkdir(path)
        if atomic_env_specific:
            path = os.path.join(path, idd)#if idd dir exists
            if not os.path.exists(path):
                os.mkdir(path)
        for f in os.listdir(path):
            nm, ex = os.path.splitext(f)
            if nm==check_name and os.path.isfile(os.path.join(path,f)):
                return True, os.path.join(path, f)  #results already exist

        return False, path
            
    def store_to_file(self, result, path, fname, file_ext=None):
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


    def generate_file_name(self, descriptor, idd, z, **args):
        sep=""
        name = sep.join([descriptor, ":", z, "_",idd,"__|"])
        for i in args:
            name = sep.join([name, i, "_", str(args[i]), "|"])

        return name

