"""Functions and Store class for storing atomic descriptions
"""
import os
import numpy as np
import pickle
import shutil
import time
import types


class Store:
    """Class for efficient storing of description of the AtomsCollection"""

    def __init__(self, store_path=None):
        """Default initiation creates a Store entitled 'store' in current directory"""
        if store_path is None:
            self.root = os.path.join(os.getcwd(), "store")
        else:
            self.root = os.path.expanduser(store_path)
        if not os.path.exists(self.root):
            os.mkdir(self.root)

    def __str__(self):
        """Returns relative path of the store location for the string representation of Store object"""
        return os.path.relpath(self.root)

    def list_store_results(self):
        """Print results currently stored in the Store"""
        #TODO
        pass

    def _generate_default_file_name(self, param_1, param_2):
        """Generate a default file name based on given parameters and the current date and time."""
        timestr = time.strftime("%Y%m%d-%H%M%S")
        return param_1 + "_" + param_2 + "_" + timestr + ".pkl"

    def _store_file(self, to_store, path):
        """Pickle and store result to the given path"""
        with open(path, 'wb') as f:
            pickle.dump(to_store, f)

    def store_description(self, result, aid, descriptor, info, **desc_args):
        """Function to store information into result store

        Parameters:
            result (pickle): computed result
            aid (str): atoms id
            descriptor (str): name of descriptor
            info (dict): dictionary holding additional information to be stored, including the parameters used in calculation
            kwargs (dict): arguments for computing descriptor, will be used to generate file names

        """
        fname = self._generate_default_file_name(aid, descriptor)
        info_fname = "info_" + fname

        path = os.path.join(self.root, "Descriptions", aid, descriptor)
        os.makedirs(path, exist_ok=True)

        full_path = os.path.join(path, fname)
        self._store_file(result, full_path)

        #put description args in info dict
        info["desc_args"] = desc_args

        #store info dict
        info_path = os.path.join(path, info_fname)
        self._store_file(info, info_path)

    def store_collection_result(self, result, info, method, collection_name, based_on, **method_args):
        """Store collection specific results generated.

        Parameters:
            result: result to be stored
            info (dict): additional info dict to be stored simultaneously and used in result retrieval
            method (str): string giving the method used in generating results
            collection_name (str): name of the collection
            based_on (tuple of type (str, dict)): tuple containting a) the name of the descriptor the method is based
            on, and b) a dictionary holding the arguments used in the calculation of the descriptor results
            **method_args: additional arguments to used in generating the method result
        """
        fname = self._generate_default_file_name(collection_name, method)
        info_fname = "info_" + fname
        path = os.path.join(self.root, "Collections", collection_name, method)
        os.makedirs(path, exist_ok=True)

        #store result
        full_path = os.path.join(path, fname)
        self._store_file(result, full_path)

        #put description and method args in info dict
        info["based_on_name"] = based_on[0]
        info["based_on_args"] = based_on[1]
        info["method_args"] = method_args

        #store info dict
        info_path = os.path.join(path, info_fname)
        self._store_file(info, info_path)

    def store_additional(self, result, store_as, info=None):
        """Function to store any arbitrary results specified by the user"""
        #TODO
        pass

    def _unpickle(self, path, filename=None):
        """Unpickle the result stored at the given file path. Can give path to directory and file name separately,
        or just pass a single path with all the information.

        Parameters:
            path (str): string holding the file path
            filename (str): name of the file (if not already included in "path"
        """
        if filename is not None:
            path = os.path.join(path, filename)
        result = None
        try:
            with open(path, 'rb') as f:
                result = pickle.load(f)
        except FileNotFoundError as e:
            raise e
        except pickle.UnpicklingError:
            raise pickle.UnpicklingError("Exception when loading file %s, consider deleting result and recomputing" % (filename))
        return result

    def _equal_args(self, args_given, args_store):
        """Used in _get_collection_result to compare two argument dictionaries, specifically, it considers two functions the same
        if they have the same name.
        """
        if args_given == args_store:
            return True
        if len(args_given) != len(args_store):
            return False

        for key in args_given.keys():
            try:
                item_given = args_given[key]
                item_store = args_store[key]
            except:
                return False

            if isinstance(item_given, types.FunctionType) and isinstance(item_store, types.FunctionType):
                if item_given.__name__ != item_store.__name__:
                    return False
            else:
                if item_given != item_store:
                    return False
        return True

    def check_exists(self, store_section, level1, level2, based_on=None, explicit=False, **kwargs):
        """ Function to check if correct file structure is in place and if a result file exists for these parameters

        Parameters:
            store_section (str): corresponds to the store section being considered, the string "Descriptions" or "Collections".
            level1 (str): In "Descriptions" section, corresponds to the aid. In the "Collections" section, corresponds to the method name.
            level2 (str): In "Descriptions" section, corresponds to the descriptor name. In the "Collections" section, corresponds to the
            collection name.
            based_on (string, dict): If based_on is None, considered an atomic description. If not, considered a collection processing result.
        """
        path = os.path.join(self.root, store_section, level1, level2)
        if os.path.exists(path):
            directory = os.fsencode(path)

            for file in os.listdir(directory):
                filename = os.fsdecode(file)
                print(level1 +"_" + level2, filename[:-20])
                if (level1 + "_" + level2) == filename[:-20]: #remove the date/time and file end
                    info = self._unpickle(path, "info_" + filename)
                    check_args = info['desc_args'] if 'desc_args' in info else info['method_args']
                    print("Based on: ",self._based_on_is_correct(based_on, info))
                    print("Args: ", self._equal_args(kwargs, check_args))
                    if self._based_on_is_correct(based_on, info) and self._equal_args(kwargs, check_args):
                        if explicit:
                            return filename
                        else:
                            return True
            else:
                return False
        else:
            return False


    def _based_on_is_correct(self, based_on, info):
        """Function used when checking if result exists, makes sure what the user expects for the 'based_on' parameter
        matches with what is in the info dictionary."""
        #if based_on is None, and not in dict, true
        #if based_on is None, and is in dict, false
        #if based_on is not None, and not in dict or wrong in dict, false
        #if based_on is not None, and info matches, true

        name = None
        args = None
        try:
            name = info['based_on_name']
            args = info['based_on_args']
        except KeyError: #not found in dict, key error
            print("a")
            pass

        if based_on is None and (name is not None or args is not None):
            print("b")
            print(f"based_on={based_on}, name={name}, args={args}")
            return False
        elif based_on is not None:
            if (name is None or args is None) or not self._equal_args(based_on[1], args) or name != based_on[0]:
                print("c")
                return False

        return True

    def get_description(self, aid, descriptor, **desc_args):
        """Function to retrieve description results from the store for the given atoms id and parameters.

        Parameters:
            aid (str): atoms id
            descriptor (str): descriptor name
            **desc_args: arguments used in generating the result
        """
        filename = self.check_exists("Descriptions", aid, descriptor, explicit=True, **desc_args)
        if filename is False:
            raise FileNotFoundError("No such results found for given parameters")
        else:
            path = os.path.join(self.root, "Descriptions", aid, descriptor)
            res = self._unpickle(path, filename)
            info = self._unpickle(path, "info_" + filename)
            return res, info

    def get_collection_result(self, method, collection_name, based_on, **method_args):
        filename = self.check_exists("Collections", collection_name, method, based_on=based_on, explicit=True, **method_args)
        if filename is False:
            raise FileNotFoundError("No such results found for given parameters")
        else:
            path = os.path.join(self.root, "Collections", collection_name, method)
            res = self._unpickle(path, filename)
            info = self._unpickle(path, "info_" + filename)
            return res, info

    def _clear_result(self, descriptor, idd, **kwargs):
        '''Function to remove a single result from the store

        Parameters:
            descriptor (str): name of descriptor
            idd (str): atoms id
            kwargs (dict): refers to whatever arguments that were used to generate the description (which correspond to the file name)
        '''
        #FIXME: Store structure changed
        # fname = self._generate_file_name(descriptor, idd, **kwargs)
        # path = os.path.join(self.root, descriptor, idd, fname)
        # if os.path.exists(path):
        #     os.remove(path)
        # directory = os.path.dirname(path)
        # if os.path.isdir(directory) and len(os.listdir(directory)) == 0:
        #     os.rmdir(directory)
        pass

    def clear(self, descriptor, idd, **kwargs):
        '''More general clear function

        Parameters:
            descriptor (str): name of descriptor
            idd (str **or** list(str)): atoms id (or list of atom id's)
            kwargs (dict): refers to whatever arguments that were used to generate the description (which correspond to the file name)

        '''
        # FIXME: Store structure changed
        # if type(idd) is list:
        #     for i in idd:
        #         self._clear_result(descriptor, i, **kwargs)
        # else:
        #     self._clear_result(descriptor, idd, **kwargs)
        # directory = os.path.join(self.root, descriptor)
        # if os.path.isdir(directory) and len(os.listdir(directory)) == 0:
        #     os.rmdir(directory)
        pass

    def clear_descriptor(self, descriptor):
        '''Function to remove all descriptions of a certain type

        Parameters:
            descriptor (str): name of descriptor

        '''
        # FIXME: Store structure changed
        # path = os.path.join(self.root, descriptor)
        # if os.path.exists(path):
        #     shutil.rmtree(path)
        pass

    def clear_all(self):
        '''Function to remove all results from the Store'''
        # FIXME: Store structure changed
        # for item in os.listdir(self.root):
        #     path = os.path.join(self.root, item)
        #     if os.path.isdir(path):
        #         shutil.rmtree(path)
        pass
