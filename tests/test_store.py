from gblearn.store import Store as rs
import os
import numpy as np
import shutil
import unittest


class TestResultStore(unittest.TestCase):

    def test_check_exists(self):
        '''Test check_existing_results function'''
        r1 = rs('./tests/test_paths/')
        e1 = r1.check_exists('desc', 'aid1', a="result")
        assert e1 == True
        e2 = r1.check_exists('desc', 'aid1', a="b")
        assert e2 == False
        e3 = r1.check_exists('desc', 'aid2', rcut=9.0, nmax=11, lmax=11)
        assert e3 == False
        e4 = r1.check_exists('desc', 'no_file', a="n")
        assert e4 == False
        e5 = r1.check_exists('no_file', 'aid1', a='result1')
        assert e5 == False

    def test_generate_file_name(self):
        '''Test generate file name function'''
        r2 = rs(".")
        f1 = r2._generate_file_name("soap", "111", rcut=9.0, nmax=11, lmax=11)
        assert f1 == "soap__111___lmax_11___nmax_11___rcut_9.0.pkl"
        f2 = r2._generate_file_name("soap", "111",  nmax=11, rcut=9, lmax=11)
        assert f2 == "soap__111___lmax_11___nmax_11___rcut_9.pkl"
        f3 = r2._generate_file_name("soap", "111",  nmax=11, rcut=9.0, lmax=11)
        assert f3 == "soap__111___lmax_11___nmax_11___rcut_9.0.pkl"
        f4 = r2._generate_file_name(
            "ler", 'U', collection="A", eps=0.025, rcut=5.0, nmax=9,
            lmax=9, metric="euclidean", n_trees=10, search_k=-1)
        assert f4 == "ler__U___collection_A___eps_0.025___lmax_9___metric_euclidean___n_trees_10___nmax_9___rcut_5.0___search_k_-1.pkl"

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

    def test_get(self):
        '''Tests to make sure get_descriptor returns expected value'''
        desc = "soap"
        aid = "aid_111"
        r4 = rs("./tests/results/")
        res = np.array([[1, 2, 3], [4, 5, 6]])
        r4.store(res, desc, aid, rcut=9, nmax=10, lmax=10)
        ret_val = r4.get(desc, aid, rcut=9, nmax=10, lmax=10)
        assert np.array_equal(ret_val, res)

        # No corresponding result (missing parameter)
        ret_val2 = r4.get(aid, desc, rcut=9, lmax=10)
        assert ret_val2 == None
        shutil.rmtree("./tests/results/")

        r3 = rs("./tests/results")
        result = "Random test result"
        desc = "test"
        aid = "12"
        r3.store(result, desc, aid, arg_a=1, arg_b=2)
        ret_val3 = r3.get(desc, aid, arg_a=1, arg_b=2)
        assert ret_val3 == result
