from pyrelate.store import Store
from pyrelate.collection import AtomsCollection
import os
import io
import sys
import numpy as np
import shutil
import unittest
# FIXME clean up and simplify unit tests


def _delete_store(store):
    shutil.rmtree(store.root)


def _test_descriptor(atoms, **kwargs):
    return 'test result'


def _initialize_collection_and_read(aids):
    '''initialize collection and read specified atoms files'''
    my_col = AtomsCollection("Test", "tests/results")
    data_path = 'tests/test_data/ni.p{0:s}.out'
    for aid in aids:
        my_col.read(data_path.format(aid),
                    28, 'lammps-dump-text', rxid=r'ni.p(?P<aid>\d+).out')
    return my_col


def _initialize_collection_and_describe(desc, aids, **kwargs):
    '''initialize collection with specified descriptor and aids and describe'''
    my_col = _initialize_collection_and_read(aids)
    for d in desc:
        my_col.describe(d, fcn=_test_descriptor, **kwargs)
        for aid in aids:
            assert my_col.get(d, aid, **kwargs) != None
    return my_col


class TestStore(unittest.TestCase):

    # test initialization
    def test_init_store(self):
        try:
            r = Store('./tests/test_paths/')
        except:
            assert False
        assert True

    # test default initialization
    def test_init_store_default(self):
        try:
            r = Store()
        except:
            assert False
        shutil.rmtree("./store")
        assert True

    # check_results function
    def test_check_exists_true(self):
        r1 = Store('./tests/test_paths/')
        e1 = r1.check_exists('desc', 'aid1', a="result")
        assert e1 == True

    def test_check_exists_false(self):
        r1 = Store('./tests/test_paths/')
        e1 = r1.check_exists('desc', 'aid1', a="b")
        assert e1 == False
        e2 = r1.check_exists('no_file', 'aid1', a='result1')
        assert e2 == False

    # test generate_file_name function
    def test_generate_file_name(self):
        r1 = Store(".")
        f1 = r1._generate_file_name("soap", "111", rcut=9.0, nmax=11, lmax=11)
        assert f1 == "soap__111___lmax_11___nmax_11___rcut_9.0.pkl"

    def test_generate_file_name_diff_order(self):
        r1 = Store(".")
        f1 = r1._generate_file_name("soap", "111",  nmax=11, rcut=9.0, lmax=11)
        assert f1 == "soap__111___lmax_11___nmax_11___rcut_9.0.pkl"

    def test_generate_file_name_with_collection(self):
        my_col = AtomsCollection("A", ".")
        f = my_col.store._generate_file_name(
            "ler", 'U', collection=my_col, eps=0.025, rcut=5.0, nmax=9,
            lmax=9, metric="euclidean", n_trees=10, search_k=-1)
        assert f == "ler__U___collection_a___eps_0.025___lmax_9___metric_euclidean___n_trees_10___nmax_9___rcut_5.0___search_k_-1.pkl"

    def test_store(self):
        '''Tests storing results as a numpy array or other type'''
        r3 = Store("./tests/results")
        result = "Random test result"
        desc = "test_desc"
        aid = "111"
        kw1 = "option_1"
        kw2 = "option_2"
        r3.store(result, desc, aid, a=kw1, b=kw2)
        fname = r3._generate_file_name(desc, aid, a=kw1, b=kw2)
        fpath = os.path.join(r3.root, desc, aid, fname)
        assert os.path.exists(fpath)
        shutil.rmtree("./tests/results/")

    def test_get_file_unpickling_error(self):
        desc = "desc"
        aid = "fakepkl"
        r = Store("./tests/test_paths/")
        output = io.StringIO()
        sys.stdout = output
        r.get(desc, aid, arg1="1")
        assert "UnpicklingError when loading file desc__fakepkl___arg1_1.pkl, consider deleting result and recomputing\n" == output.getvalue()

    def test_get_file_numpy_array(self):
        '''Tests to make sure get_descriptor returns expected value'''
        desc = "soap"
        aid = "aid_111"
        r4 = Store("./tests/results/")
        res = np.array([[1, 2, 3], [4, 5, 6]])
        r4.store(res, desc, aid, rcut=9, nmax=10, lmax=10)
        ret_val = r4._get_file(desc, aid, rcut=9, nmax=10, lmax=10)
        assert np.array_equal(ret_val, res)
        shutil.rmtree("./tests/results/")

        # No corresponding result (missing parameter)
    def test_get_file_missing_param(self):
        desc = "soap"
        aid = "aid_111"
        r4 = Store("./tests/results/")
        res = np.array([[1, 2, 3], [4, 5, 6]])
        r4.store(res, desc, aid, rcut=9, nmax=10, lmax=10)
        ret_val2 = r4._get_file(aid, desc, rcut=9, lmax=10)
        assert ret_val2 == None
        shutil.rmtree("./tests/results/")

    def test_get_file_string(self):
        r3 = Store("./tests/results")
        result = "Random test result"
        desc = "test"
        aid = "12"
        r3.store(result, desc, aid, arg_a=1, arg_b=2)
        ret_val3 = r3._get_file(desc, aid, arg_a=1, arg_b=2)
        assert ret_val3 == result
        shutil.rmtree("./tests/results/")

    def test_get(self):
        r3 = Store("./tests/test_paths")
        res = r3.get("soap", '455', rcut=5.0, nmax=9, lmax=9)
        assert type(res) is np.ndarray

    def test_get_with_list(self):
        r3 = Store("./tests/test_paths")
        res = r3.get("soap", ['455'], rcut=5.0, nmax=9, lmax=9)
        assert type(res) == dict
        assert type(res['455']) is not type(None)
        assert len(res) == 1

    def test_clear_result(self):
        my_col = _initialize_collection_and_describe(
            ['test'], ['454', '455'], arg1=1)
        my_col.store._clear_result('test', '454', arg1=1)
        assert my_col.store.check_exists('test', '454', arg1=1) == False
        assert my_col.store.check_exists('test', '455', arg1=1) == True
        assert os.path.exists(os.path.join(
            my_col.store.root, 'test', '454')) == False
        assert os.path.exists(my_col.store.root) == True

    def test_clear_specific_results_for_collection(self):
        my_col = _initialize_collection_and_describe(
            ['test'], ['454', '455'], arg1=1)
        my_col.store.clear('test', ['454', '455'], arg1=1)
        assert my_col.store.check_exists('test', '454', arg1=1) == False
        assert my_col.store.check_exists('test', '455', arg1=1) == False
        assert os.path.exists(os.path.join(
            my_col.store.root, 'test', '454')) == False
        assert os.path.exists(os.path.join(my_col.store.root, 'test')) == False
        assert os.path.exists(my_col.store.root) == True

    def test_clear_descriptor(self):
        my_col = _initialize_collection_and_describe(
            ['test', 'test2'], ['454'], arg1=1)
        my_col.store.clear_descriptor('test')
        assert my_col.store.check_exists('test', '454', arg1=1) == False
        assert my_col.store.check_exists('test2', '454', arg1=1) == True
        assert os.path.exists(os.path.join(my_col.store.root, 'test')) == False

    def test_clear_all(self):
        my_col = _initialize_collection_and_describe(
            ['test', 'test2'], ['454', '455'], arg1=1)
        my_col.store.clear_all()
        assert os.path.exists(os.path.join(
            my_col.store.root, 'test', '454')) == False
        assert os.path.exists(os.path.join(my_col.store.root, 'test')) == False
        assert os.path.exists(os.path.join(
            my_col.store.root, 'test2')) == False
        assert os.path.exists(my_col.store.root) == True
