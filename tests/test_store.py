from pyrelate.store import Store as rs
from pyrelate.collection import AtomsCollection
import os
import numpy as np
import shutil
import unittest


class TestStore(unittest.TestCase):

    #test initialization
    def test_init_store(self):
        try:
            r = rs('./tests/test_paths/')
        except:
            assert False
        assert True

    #test default initialization
    def test_init_store_default(self):
        try:
            r = rs()
        except:
            assert False
        shutil.rmtree("./store")
        assert True

    #check_results function
    def test_check_exists_true(self):
        r1 = rs('./tests/test_paths/')
        e1 = r1.check_exists('desc', 'aid1', a="result")
        assert e1 == True

    def test_check_exists_false(self):
        r1 = rs('./tests/test_paths/')
        e1 = r1.check_exists('desc', 'aid1', a="b")
        assert e1 == False
        e2 = r1.check_exists('no_file', 'aid1', a='result1')
        assert e2 == False

    #test generate_file_name function
    def test_generate_file_name(self):
        r1 = rs(".")
        f1 = r1._generate_file_name("soap", "111", rcut=9.0, nmax=11, lmax=11)
        assert f1 == "soap__111___lmax_11___nmax_11___rcut_9.0.pkl"

    def test_generate_file_name_diff_order(self):
        r1 = rs(".")
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
        r3 = rs("./tests/results")
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

    def test_get_file_numpy_array(self):
        '''Tests to make sure get_descriptor returns expected value'''
        desc = "soap"
        aid = "aid_111"
        r4 = rs("./tests/results/")
        res = np.array([[1, 2, 3], [4, 5, 6]])
        r4.store(res, desc, aid, rcut=9, nmax=10, lmax=10)
        ret_val = r4._get_file(desc, aid, rcut=9, nmax=10, lmax=10)
        assert np.array_equal(ret_val, res)
        shutil.rmtree("./tests/results/")

        # No corresponding result (missing parameter)
    def test_get_file_missing_param(self):
        desc = "soap"
        aid = "aid_111"
        r4 = rs("./tests/results/")
        res = np.array([[1, 2, 3], [4, 5, 6]])
        r4.store(res, desc, aid, rcut=9, nmax=10, lmax=10)
        ret_val2 = r4._get_file(aid, desc, rcut=9, lmax=10)
        assert ret_val2 == None
        shutil.rmtree("./tests/results/")

    def test_get_file_string(self):
        r3 = rs("./tests/results")
        result = "Random test result"
        desc = "test"
        aid = "12"
        r3.store(result, desc, aid, arg_a=1, arg_b=2)
        ret_val3 = r3._get_file(desc, aid, arg_a=1, arg_b=2)
        assert ret_val3 == result
        shutil.rmtree("./tests/results/")

    def test_get(self):
        r3 = rs("./tests/test_paths")
        res = r3.get("soap", '455', rcut=5.0, nmax=9, lmax=9)
        assert type(res) is np.ndarray

    def test_get_with_list(self):
        r3 = rs("./tests/test_paths")
        res = r3.get("soap", ['455'], rcut=5.0, nmax=9, lmax=9)
        assert type(res) == dict
        assert type(res['455']) is not type(None)
        assert len(res) == 1
