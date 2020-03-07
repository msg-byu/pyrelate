from gblearn.store import ResultStore as rs
import os
import numpy as np
import shutil
import unittest


class TestResultStore(unittest.TestCase):

    def test_check_existing_results(self):
        '''Test check_existing_results function
        '''
        r1 = rs('./tests/test_paths/')
        e1 = r1.check_existing_results('desc', 'aid1', 'result1')
        assert e1 == True
        e2 = r1.check_existing_results('desc', 'aid1', 'nonex')
        assert e2 == False
        e3 = r1.check_existing_results('desc', 'aid2', 'result1')
        assert e3 == False
        e4 = r1.check_existing_results('desc', 'no_file', 'result1')
        assert e4 == False
        e5 = r1.check_existing_results('no_file', 'aid1', 'result1')
        assert e5 == False

    def test_generate_file_name(self):
        '''Test generate file name function'''
        r2 = rs(".")
        f1 = r2.generate_file_name("soap", "111", rcut=9, nmax=11, lmax=11)
        assert f1 == "soap__111___lmax_11___nmax_11___rcut_9.0"
        f2 = r2.generate_file_name("soap", "111",  nmax=11, rcut=9, lmax=11)
        assert f2 == "soap__111___lmax_11___nmax_11___rcut_9.0"
        f3 = r2.generate_file_name("soap", "111",  nmax=11, rcut=9.0, lmax=11)
        assert f3 == "soap__111___lmax_11___nmax_11___rcut_9.0"
        # test should throw value error if rcut (for SOAP descriptor) is not an int or float
        self.assertRaises(ValueError, r2.generate_file_name, "soap", "111",
                          nmax=11, rcut="a", lmax=11)

    def test_store_descriptor(self):
        '''Tests storing results as a numpy array or other type'''
        # test storing random result type
        r3 = rs("./tests/results")
        result = "Random test result"
        desc = "test_desc"
        aid = "aid_111"
        fname = "filename"
        r3.store_descriptor(result, desc, aid, fname)
        fname = fname + ".dat"
        fpath = os.path.join(r3.root, desc, aid, fname)
        assert os.path.exists(fpath)
        shutil.rmtree("./tests/results/")

        # test storing numpy ndarray
        r4 = rs("./tests/results")
        result2 = np.array([[1, 2, 3], [4, 5, 6]])
        r4.store_descriptor(result2, desc, aid, fname)
        fname = fname + ".npy"
        fpath = os.path.join(r4.root, desc, aid, fname)
        assert os.path.exists(fpath)
        shutil.rmtree("./tests/results/")

    def test_get_descriptor(self):
        '''Tests to make sure get_descriptor returns expected value'''
        desc = "soap"
        aid = "aid_111"
        fname = "soap__aid_111___lmax_10___nmax_10___rcut_9.0"
        r4 = rs("./tests/results/")
        res = np.array([[1, 2, 3], [4, 5, 6]])
        r4.store_descriptor(res, desc, aid, fname)
        ret_val = r4.get_descriptor(aid, desc, rcut=9, nmax=10, lmax=10)
        assert np.array_equal(ret_val, res)
        # missing parameters
        ret_val2 = r4.get_descriptor(aid, desc, rcut=9, nmax=10)
        assert ret_val2 == -1
        shutil.rmtree("./tests/results/")

        r3 = rs("./tests/results")
        result = "Random test result"
        desc = "test"
        aid = "12"
        fname = r3.generate_file_name(desc, aid, arg_a=1, arg_b=2)
        r3.store_descriptor(result, desc, aid, fname)
        ret_val3 = r3.get_descriptor(aid, desc, arg_a=1, arg_b=2)
        assert ret_val3 == result
