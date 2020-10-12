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

    def test_check_exists_true(self):
        '''Test check_exists, result does exist'''
        store = Store('./tests/test_paths/')
        exists = store.check_exists('desc', 'aid1', a="result")
        assert exists == True

    def test_check_exists_wrong_kwargs(self):
        '''Test check_exists, result with given params doesn't exist'''
        store = Store('./tests/test_paths/')
        exists = store.check_exists('desc', 'aid1', a="b")
        assert exists == False

    def test_check_exists_wrong_descriptor(self):
        '''Test check_exists, result with given descriptor doesn't exist'''
        store = Store('./tests/test_paths/')
        exists = store.check_exists('wrong_desc', 'aid1', a='result1')
        assert exists == False

    def test_check_exists_wrong_aid(self):
        '''Test check_exists, result for given aid doesn't exist'''
        store = Store('./tests/test_paths/')
        exists = store.check_exists('desc', 'wrong_aid', a='result1')
        assert exists == False

    def test_generate_file_name(self):
        '''Test _generate_file_name'''
        store = Store("./tests/results")
        filename = store._generate_file_name(
            "soap", "111", rcut=9.0, nmax=11, lmax=11)
        assert filename == "soap__111___lmax_11___nmax_11___rcut_9.0.pkl"
        _delete_store(store)

    def test_generate_file_name_diff_order(self):
        '''Test _generate_file_name, check same filename for parameters in different order'''
        store = Store("./tests/results")
        filename = store._generate_file_name(
            "soap", "111",  nmax=11, rcut=9.0, lmax=11)
        assert filename == "soap__111___lmax_11___nmax_11___rcut_9.0.pkl"
        _delete_store(store)

    def test_generate_file_name_with_collection(self):
        '''Test _generate_file_name with collection as one parameter (make sure __str__ method called)'''
        my_col = AtomsCollection("A", "./tests/results")
        filename = my_col.store._generate_file_name(
            "ler", 'U', collection=my_col, eps=0.025, rcut=5.0, nmax=9, lmax=9, metric="euclidean", n_trees=10, search_k=-1)
        assert filename == "ler__U___collection_a___eps_0.025___lmax_9___metric_euclidean___n_trees_10___nmax_9___rcut_5.0___search_k_-1.pkl"
        _delete_store(my_col.store)

    def test_generate_file_name_with_function(self):
        '''Test _generate_file_name with function as one of the inputs'''
        my_col = AtomsCollection("A", "./tests/results")
        filename = my_col.store._generate_file_name(
            "desc", 'aid', fcn_required = _test_descriptor)
        assert filename == "desc__aid___fcn_required__test_descriptor.pkl"
        _delete_store(my_col.store)


    def test_store(self):
        '''Tests storing results as a pickle'''
        store = Store("./tests/results")
        result = "Random test result"
        desc = "test_desc"
        aid = "111"
        kw1 = "option_1"
        kw2 = "option_2"
        store.store(result, desc, aid, a=kw1, b=kw2)
        fname = store._generate_file_name(desc, aid, a=kw1, b=kw2)
        fpath = os.path.join(store.root, desc, aid, fname)
        assert os.path.exists(fpath)
        _delete_store(store)

    def test_get_file_unpickling_error(self):
        '''Test get, unpickling error'''
        desc = "desc"
        aid = "fakepkl"
        store = Store("./tests/test_paths/")
        output = io.StringIO()
        sys.stdout = output
        store.get(desc, aid, arg1="1")
        assert "UnpicklingError when loading file desc__fakepkl___arg1_1.pkl, consider deleting result and recomputing\n" == output.getvalue()

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
