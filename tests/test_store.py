from pyrelate.store import Store
from pyrelate.collection import AtomsCollection
import os
import io
import sys
import numpy as np
import shutil
import unittest
#TODO finish updating unit tests

def _delete_store(store):
    shutil.rmtree(store.root)


def _test_descriptor(atoms, **kwargs):
    return 'test result'


def _initialize_collection_and_read(aids):
    '''Initialize collection and read specified atoms files'''
    my_col = AtomsCollection("Test", "tests/results")
    data_path = 'tests/test_data/ni.p{0:s}.out'
    for aid in aids:
        my_col.read(data_path.format(aid), 28, 'lammps-dump-text',
                    rxid=r'ni.p(?P<aid>\d+).out')
    return my_col


def _initialize_collection_and_describe(desc, aids, **kwargs):
    '''Initialize collection with specified descriptor and aids and describe'''
    my_col = _initialize_collection_and_read(aids)
    for d in desc:
        my_col.describe(d, fcn=_test_descriptor, **kwargs)
        for aid in aids:
            assert my_col.get(d, aid, **kwargs) != None
    return my_col


class TestStore(unittest.TestCase):

    def test_init_store(self):
        '''Test initialization of store'''
        try:
            store = Store('./tests/results')
        except:
            assert False
        assert True
        _delete_store(store)

    def test_init_store_default(self):
        '''Test initializing store default'''
        try:
            store = Store()
        except:
            assert False
        assert True
        _delete_store(store)

    def test_init_store_expanduser(self):
        '''Test initializing store default'''
        try:
            store = Store("~/test_store")
        except:
            assert False
        assert True
        _delete_store(store)

    #Functions to test:
    ## check_exists (when done writing)
    ## store_additional (when written)
    ## _equal_args
    ## _get_description
    ## _get_collection_result

    #Done
    ## _generate_default_file_name
    ## _store_file
    ## store_description
    ## store_collection_result
    ## _unpickle

    def test_based_on_is_correct_1(self):
        # if based_on is None, and not in dict, true
        store = Store("./tests/results")
        based_on = None
        info = {"other": 1}
        assert store._based_on_is_correct(based_on, info)
        _delete_store(store)

    def test_based_on_is_correct_2(self):
        # if based_on is None, and is in dict, false
        store = Store("./tests/results")
        based_on = None
        info = {"other": 1, "based_on_name":"name", "based_on_args":{"a":1,"b":2}}
        assert not store._based_on_is_correct(based_on, info)
        _delete_store(store)

    def test_based_on_is_correct_3(self):
        # if based_on is not None, and not in dict, false
        store = Store("./tests/results")
        desc_args = {"a":1,"b":2}
        based_on = ("desc", desc_args)
        info = {"other": 1}
        assert not store._based_on_is_correct(based_on, info)
        _delete_store(store)

    def test_based_on_is_correct_4(self):
        # if based_on is not None, and wrong based_on_args, false
        store = Store("./tests/results")
        desc_args = {"a":1,"b":2}
        name = "desc"
        based_on = (name, desc_args)
        info = {"other": 1, "based_on_name":name, "based_on_args":{"a":1,"b":22}}
        assert not store._based_on_is_correct(based_on, info)
        _delete_store(store)

    def test_based_on_is_correct_5(self):
        # if based_on is not None, and info matches, true
        store = Store("./tests/results")
        desc_args = {"a": 1, "b": 2}
        name = "desc"
        based_on = (name, desc_args)
        info = {"other": 1, "based_on_name": name, "based_on_args": desc_args}
        assert store._based_on_is_correct(based_on, info)
        _delete_store(store)

    def test_check_exists_description_true(self):
        store = Store("./tests/results")
        result = "Random test result"
        desc = "test_desc"
        aid = "111"
        kw1 = "option_1"
        kw2 = "option_2"
        info = {}
        store.store_description(result, aid, desc, info, a=kw1, b=kw2)

        assert store.check_exists("Descriptions", aid, desc, a=kw1, b=kw2)
        assert type(store.check_exists("Descriptions", aid, desc, a=kw1, b=kw2, explicit=True)) is str
        _delete_store(store)

    def test_check_exists_description_false(self):
        store = Store("./tests/results")
        desc = "test_desc"
        aid = "111"
        kw1 = "option_1"
        kw2 = "option_2"

        assert not store.check_exists("Descriptions", aid, desc, a=kw1, b=kw2)
        _delete_store(store)

    def test_check_exists_collection_true(self):
        store = Store("./tests/results")
        result = "Random test result"
        info = {
            "additional_info": [1, 2, 3, 4, 5]
        }
        desc = "test_desc"
        method = "test_method"
        name = "my_collection"
        desc_args = {
            "kw1": "option_1",
            "kw2": "option_2"
        }
        method_args = {
            "eps": 1,
            "num": 50
        }

        store.store_collection_result(result, info, method, name, (desc, desc_args), **method_args)

        assert store.check_exists("Collections", name, method, based_on=(desc, desc_args), **method_args)
        _delete_store(store)

    def test_check_exists_collection_false(self):
        store = Store("./tests/results")
        desc = "test_desc"
        method = "test_method"
        name = "my_collection"
        desc_args = {
            "kw1": "option_1",
            "kw2": "option_2"
        }
        method_args = {
            "eps": 1,
            "num": 50
        }

        assert not store.check_exists("Collections", name, method, based_on=(desc, desc_args), **method_args)
        _delete_store(store)

    def test_generate_default_file_name(self):
        '''Test _generate_default_file_name'''
        store = Store("./tests/results")
        aid = "455"
        desc = "soap"
        filename = store._generate_default_file_name(aid, desc)
        assert filename[:-20] == aid + "_" + desc
        assert filename[-4:] == ".pkl"
        _delete_store(store)

    def test_store_file(self):
        test = "thing"
        store = Store("./tests/results")
        path = os.path.join(store.root, "thing.pkl")
        store._store_file(test, path)
        assert os.path.exists(path)
        _delete_store(store)

    def test_store_description(self):
        #called in describe()
        #result, descriptor, aid, argmuments
        store = Store("./tests/results")
        result = "Random test result"
        desc = "test_desc"
        aid = "111"
        kw1 = "option_1"
        kw2 = "option_2"
        info = {}
        store.store_description(result, aid, desc, info, a=kw1, b=kw2)

        fpath = os.path.join(store.root, "Descriptions", aid, desc)
        info_fpath = os.path.join(store.root, "Descriptions", aid, desc)
        assert os.path.exists(fpath)

        directory = os.fsencode(fpath)
        for file in os.listdir(directory):
            filename = os.fsdecode(file)
            if filename[:-20] == aid + "_" + desc:
                fullpath = os.path.join(fpath, filename)
                assert os.path.exists(fullpath)
                break
        else:
            assert False, "No correct file found"
        _delete_store(store)

    def test_store_description_with_info(self):
        # called in describe()
        # result, descriptor, aid, argmuments
        store = Store("./tests/results")
        result = "Random test result"
        desc = "test_desc"
        aid = "111"
        kw1 = "option_1"
        kw2 = "option_2"
        info = {"num":47, "important_info":12 }
        store.store_description(result, aid, desc, info=info, a=kw1, b=kw2)

        fpath = os.path.join(store.root, "Descriptions", aid, desc)
        assert os.path.exists(fpath)

        directory = os.fsencode(fpath)
        for file in os.listdir(directory):
            filename = os.fsdecode(file)
            if filename[:-20] == aid + "_" + desc:
                fullpath = os.path.join(fpath, filename)
                assert os.path.exists(fullpath)
                info_fullpath = os.path.join(fpath, "info_" + filename)
                assert os.path.exists(info_fullpath)
                fetched_info = store._unpickle(info_fullpath)
                assert fetched_info['num'] == info['num']
                assert fetched_info['important_info'] == info['important_info']
                assert fetched_info['desc_args'] == {"a": kw1, "b": kw2}
                break
        else:
            assert False, "No correct file found"
        _delete_store(store)


    def test_store_collection_description(self):
        #called in the process() method
        #result, info, collection name, arguments, descriptor_args
        store = Store("./tests/results")
        result = "Random test result"
        info = {
            "additional_info": [1,2,3,4,5]
        }
        desc = "test_desc"
        method = "test_method"
        name = "my_collection"
        desc_args = {
            "kw1": "option_1",
            "kw2": "option_2"
        }
        method_args = {
            "eps": 1,
            "num":50
        }

        store.store_collection_result(result, info, method, name, (desc, desc_args), **method_args)

        fpath = os.path.join(store.root, "Collections", name, method)
        assert os.path.exists(fpath)

        directory = os.fsencode(fpath)
        for file in os.listdir(directory):
            filename = os.fsdecode(file)
            print(filename[:-20], name + "_" + method)
            if filename[:-20] == name + "_" + method:
                fullpath = os.path.join(fpath, filename)
                assert os.path.exists(fullpath)
                info_fullpath = os.path.join(fpath, "info_" + filename)
                assert os.path.exists(info_fullpath)
                fetched_info = store._unpickle(info_fullpath)
                assert fetched_info['method_args'] == method_args
                assert fetched_info['based_on_name'] == desc
                assert fetched_info['based_on_args'] == desc_args
                assert fetched_info['additional_info'] == info['additional_info']
                break
        else:
            assert False, "No correct file found"
        _delete_store(store)

    def test_store_additional(self):
        pass

    def test_unpickle_path_and_fname(self):
        store = Store("./tests/results")
        test = "thing"
        fname = "thing.pkl"
        path = os.path.join(store.root, fname)
        store._store_file(test, path)

        fetched = store._unpickle(store.root, fname)
        assert fetched == test
        _delete_store(store)

    def test_unpickle_path(self):
        store = Store("./tests/results")
        test = "thing"
        fname = "thing.pkl"
        path = os.path.join(store.root, fname)
        store._store_file(test, path)

        fetched = store._unpickle(path)
        assert fetched == test
        _delete_store(store)

    def test_unpickle_does_not_exist(self):
        store = Store("./tests/results")
        try:
            store._unpickle("fakepath", "fake_fname")
        except FileNotFoundError as e:
            assert True
        else:
            assert False, "Expected error not thrown"

    def test_unpickle_unpickling_error(self):
        import pickle
        store = Store("./tests/test_paths/")
        fname = "fakepkl.pkl"
        try:
            store._unpickle(store.root, fname)
        except pickle.UnpicklingError as e:
            assert True
        else:
            assert False, "Expected error not thrown"

    """
    def test_get_file_unpickling_error(self):
        '''Test get, unpickling error'''
        desc = "desc"
        aid = "fakepkl"
        store = Store("./tests/test_paths/")
        output = io.StringIO()
        sys.stdout = output
        store.get(desc, aid, arg1="1")
        assert "UnpicklingError when loading file fakepkl.pkl, consider deleting result and recomputing\n" == output.getvalue()

    def test_get_file_numpy_array(self):
        '''Tests _get_file, make sure get_descriptor returns expected value'''
        desc = "soap"
        aid = "aid_111"
        store = Store("./tests/results/")
        res = np.array([[1, 2, 3], [4, 5, 6]])
        store.store(res, desc, aid, rcut=9, nmax=10, lmax=10)
        ret_val = store._get_file(desc, aid, rcut=9, nmax=10, lmax=10)
        assert np.array_equal(ret_val, res)
        _delete_store(store)
        # shutil.rmtree("./tests/results/")

    def test_get_file_missing_param(self):
        '''Test _get_file, missing parameter return None'''
        desc = "soap"
        aid = "aid_111"
        store = Store("./tests/results/")
        res = np.array([[1, 2, 3], [4, 5, 6]])
        store.store(res, desc, aid, rcut=9, nmax=10, lmax=10)
        ret_val = store._get_file(aid, desc, rcut=9, lmax=10)
        assert ret_val == None
        _delete_store(store)
        # shutil.rmtree("./tests/results/")

    def test_get_file_string(self):
        '''Test _get_file, result is string'''
        store = Store("./tests/results")
        result = "Random test result"
        desc = "test"
        aid = "12"
        store.store(result, desc, aid, arg_a=1, arg_b=2)
        ret_val = store._get_file(desc, aid, arg_a=1, arg_b=2)
        assert ret_val == result
        _delete_store(store)
        # shutil.rmtree("./tests/results/")

    def test_get(self):
        '''Test get function'''
        store = Store("./tests/test_paths")
        #result is stored in test_paths
        res = store.get("soap", '455', rcut=5.0, nmax=9, lmax=9)
        assert type(res) is np.ndarray

    def test_get_with_list(self):
        '''Test get, pass in list of aids to get results for'''
        store = Store("./tests/test_paths")
        res = store.get("soap", ['455'], rcut=5.0, nmax=9, lmax=9)
        assert type(res) == dict
        assert type(res['455']) is not type(None)
        assert len(res) == 1

    def test_clear_specific_result(self):
        '''Test clear, specific result'''
        my_col = _initialize_collection_and_describe(
            ['test'], ['454', '455'], arg1=1)
        my_col.store._clear_result('test', '454', arg1=1)
        assert my_col.store.check_exists('test', '454', arg1=1) == False
        assert my_col.store.check_exists('test', '455', arg1=1) == True
        assert os.path.exists(os.path.join(
            my_col.store.root, 'test', '454')) == False
        assert os.path.exists(my_col.store.root) == True
        _delete_store(my_col.store)

    def test_clear_specific_results_list_of_aids(self):
        '''Test clear, clear all results for given discriptor and parameters for all aids in list'''
        aids = ['454', '455']
        my_col = _initialize_collection_and_describe(
            ['test'], aids, arg1=1)
        my_col.store.clear('test', aids, arg1=1)
        assert my_col.store.check_exists('test', '454', arg1=1) == False
        assert my_col.store.check_exists('test', '455', arg1=1) == False
        assert os.path.exists(os.path.join(
            my_col.store.root, 'test', '454')) == False
        assert os.path.exists(os.path.join(my_col.store.root, 'test')) == False
        assert os.path.exists(my_col.store.root) == True
        _delete_store(my_col.store)

    def test_clear_descriptor(self):
        '''Test clear, clear all results for given descriptor'''
        my_col = _initialize_collection_and_describe(
            ['test', 'test2'], ['454'], arg1=1)
        my_col.store.clear_descriptor('test')
        assert my_col.store.check_exists('test', '454', arg1=1) == False
        assert my_col.store.check_exists('test2', '454', arg1=1) == True
        assert os.path.exists(os.path.join(my_col.store.root, 'test')) == False

    def test_clear_all(self):
        '''Test clear, clear all'''
        my_col = _initialize_collection_and_describe(
            ['test', 'test2'], ['454', '455'], arg1=1)
        my_col.store.clear_all()
        assert os.path.exists(os.path.join(
            my_col.store.root, 'test', '454')) == False
        assert os.path.exists(os.path.join(my_col.store.root, 'test')) == False
        assert os.path.exists(os.path.join(
            my_col.store.root, 'test2')) == False
        assert os.path.exists(my_col.store.root) == True
"""