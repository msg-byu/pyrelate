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

    def check_exists(self, descriptor, idd, **kwargs):
        """ Function to check if correct file structure is in place and if a result file exists for these parameters

        Parameters:
            descriptor (str): name of descriptor.
            idd (str): id of item to check for.
            kwargs (dict): Parameters associated with the description function.

        Returns:
            bool: True if file results already exists, false if they do not
        """
        #TODO: adjust for new file structure
        # fname = self._generate_file_name(descriptor, idd, **kwargs)
        # path = os.path.join(self.root, descriptor, idd, fname)
        # return os.path.exists(path)
        pass

    def _generate_default_file_name(self, descriptor):
        timestr = time.strftime("%Y%m%d-%H%M%S")
        return descriptor + "_" + timestr + ".pkl"

    def _store_file(self, to_store, path):
        with open(path, 'wb') as f:
            pickle.dump(to_store, f)

    def store_description(self, result, idd, descriptor, info=None, **desc_args):
        """
        Function to store information into result store

        Parameters:
            result (pickle): computed result
            descriptor (str): name of descriptor
            idd (str): atoms id
            kwargs (dict): arguments for computing descriptor, will be used to generate file names

        """
        fname = self._generate_default_file_name(descriptor)
        info_fname = "info_" + fname

        path = os.path.join(self.root, "Descriptions", idd, descriptor)
        os.makedirs(path, exist_ok=True)

        full_path = os.path.join(path, fname)
        self._store_file(result, full_path)

        #put description and method args in info dict
        if info is None:
            info = {}
        info["desc_args"] = desc_args

        #store info dict
        info_path = os.path.join(path, info_fname)
        self._store_file(info, info_path)

    def store_collection_result(self, result, info, method, collection_name, based_on, **method_args):
        fname = self._generate_default_file_name(method)
        info_fname = "info_" + fname
        path = os.path.join(self.root, "Collections", collection_name, method)
        os.makedirs(path, exist_ok=True)

        #store result
        full_path = os.path.join(path, fname)
        self._store_file(result, full_path)

        #put description and method args in info dict
        info["desc_name"] = based_on[0]
        info["desc_args"] = based_on[1]
        info["method_args"] = method_args

        #store info dict
        info_path = os.path.join(path, info_fname)
        self._store_file(info, info_path)

    def store_additional(self, result, store_as, info=None):
        pass

    def _unpickle(self, path, filename=None):
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
        """Used in _get_collection_result"""
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

    def get_description(self, idd, descriptor, **desc_args):
        path = os.path.join(self.root, "Descriptions", idd, descriptor)
        if os.path.exists(path):
            directory = os.fsencode(path)

            for file in os.listdir(directory):
                filename = os.fsdecode(file)
                if descriptor == filename[:len(descriptor)]:
                    info = self._unpickle(path, "info_" + filename)
                    if self._equal_args(desc_args, info['desc_args']):
                        res = self._unpickle(path, filename)
                        return res, info
            else:
                raise FileNotFoundError("No such results found for given parameters")
        else:
            raise FileNotFoundError(f"Directory {path} does not exist.")

    def get_collection_result(self, method, collection_name, based_on, **method_args):
        path = os.path.join(self.root, "Collections", collection_name, method)
        if os.path.exists(path):
            directory = os.fsencode(path)

            for file in os.listdir(directory):
                filename = os.fsdecode(file)
                if method == filename[:len(method)]:
                    info = self._unpickle(path, "info_" + filename)
                    if (based_on[0] == info['based_on_name']) and self._equal_args(based_on[1], info['based_on_args']) and self._equal_args(method_args, info['method_args']):
                        res = self._unpickle(path, filename)
                        return res, info
            else:
                raise FileNotFoundError("No such results found for given parameters")
        else:
            raise FileNotFoundError(f"Directory {path} does not exist.")

    def _clear_result(self, descriptor, idd, **kwargs):
        '''Function to remove a single result from the store

        Parameters:
            descriptor (str): name of descriptor
            idd (str): atoms id
            kwargs (dict): refers to whatever arguments that were used to generate the description (which correspond to the file name)
        '''
        #FIXME: Store structure changed
        fname = self._generate_file_name(descriptor, idd, **kwargs)
        path = os.path.join(self.root, descriptor, idd, fname)
        if os.path.exists(path):
            os.remove(path)
        directory = os.path.dirname(path)
        if os.path.isdir(directory) and len(os.listdir(directory)) == 0:
            os.rmdir(directory)

    def clear(self, descriptor, idd, **kwargs):
        '''More general clear function

        Parameters:
            descriptor (str): name of descriptor
            idd (str **or** list(str)): atoms id (or list of atom id's)
            kwargs (dict): refers to whatever arguments that were used to generate the description (which correspond to the file name)

        '''
        # FIXME: Store structure changed
        if type(idd) is list:
            for i in idd:
                self._clear_result(descriptor, i, **kwargs)
        else:
            self._clear_result(descriptor, idd, **kwargs)
        directory = os.path.join(self.root, descriptor)
        if os.path.isdir(directory) and len(os.listdir(directory)) == 0:
            os.rmdir(directory)

    def clear_descriptor(self, descriptor):
        '''Function to remove all descriptions of a certain type

        Parameters:
            descriptor (str): name of descriptor

        '''
        # FIXME: Store structure changed
        path = os.path.join(self.root, descriptor)
        if os.path.exists(path):
            shutil.rmtree(path)

    def clear_all(self):
        '''Function to remove all results from the Store'''
        # FIXME: Store structure changed
        for item in os.listdir(self.root):
            path = os.path.join(self.root, item)
            if os.path.isdir(path):
                shutil.rmtree(path)
