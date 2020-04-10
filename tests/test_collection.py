"""Tests collection.py
"""
from gblearn.collection import AtomsCollection as col
import os
import numpy as np
import sys
import shutil
import io
import unittest


class TestCollection(unittest.TestCase):
    def test_read_aid(self):
        '''
        Test get_aid function with different inputs (regex, prefix, filename),
            including incorrect Regex pattern
        '''
        t1 = col("t1", "./tests/store")
        rxid = r'ni.p(?P<gbid>\d+).out'
        import re
        c_rxid = re.compile(rxid)
        fn1 = "../homer/ni.p454.out"
        a1 = t1._read_aid(fn1, c_rxid)
        assert a1 == "454"

        prefix = "Pre"
        a2 = t1._read_aid(fn1, c_rxid, prefix)
        assert a2 == "pre_454"

        prefix = ""
        a3 = t1._read_aid(fn1, c_rxid, prefix)
        assert a3 == "_454"

        a4 = t1._read_aid(fn1, None)
        assert a4 == "ni.p454.out"

        prefix = "Test"
        a5 = t1._read_aid(fn1, None, prefix)
        assert a5 == "test_ni.p454.out"

        output2 = io.StringIO()
        sys.stdout = output2
        import re
        c_rxid3 = re.compile(r'ni.p(P<gbid>\d+).out')
        a6 = t1._read_aid(fn1, c_rxid3, prefix)
        assert "Regex found no pattern. Resolving to filename as aid.\n" == output2.getvalue()
        assert a6 == "test_ni.p454.out"

        shutil.rmtree("./tests/store")

    def test_read(self):
        '''
        Test read() function with various input types (list of file paths, directories, files, etc)
        '''

        t1 = col("Test_1", "./tests/store")
        # list of input files
        t1.read(["./tests/test_data/ni.p454.out", "./tests/test_data/ni.p453.out"], 28,
                "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out',
                prefix="TEST")
        assert 2 == len(t1)
        assert "test_454" == list(t1)[0]
        # read single file
        t1.read("./tests/test_data/ni.p455.out", 28,
                "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out',
                prefix="TEST")
        assert 3 == len(t1)
        # read directory
        t1.read("./tests/test_data/sub1/", 28,
                "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out',
                prefix="TEST")
        assert 5 == len(t1)
        # empty directory / directory and file
        t1.read(["./tests/test_data/ni.p456.out", "./tests/test_data/empty"], 28,
                "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out',
                prefix="TEST")
        assert 6 == len(t1)
        # empty list
        t1.read([], 28,
                "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out',
                prefix="TEST")
        assert 6 == len(t1)
        # will not read previously read file
        t1.read("./tests/test_data/ni.p456.out", 28,
                "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out',
                prefix="TEST")
        assert 6 == len(t1)
        # non-existent directory
        output = io.StringIO()
        sys.stdout = output
        t1.read("definitely_wrong", 28,
                "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out',
                prefix="TEST")
        assert "Invalid file path, definitely_wrong was not read.\n" == output.getvalue()
        # no file type
        self.assertRaises(StopIteration, col.read, t1, "./tests/test_data/ni.p457.out",
                          28, rxid=r'ni.p(?P<gbid>\d+).out', prefix="TEST")
        # missing parameter
        self.assertRaises(StopIteration, col.read, t1, "./tests/test_data/ni.p456.out",
                          "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out', prefix="TEST")
        # TODO test for a list of atomic numbers
        t2 = col("Test_2", "./tests/store")
        t2.read(["./tests/test_data/ni.p454.out", "./tests/test_data/ni.p453.out"], [28, 28],
                "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out',
                prefix="TEST")
        assert 2 == len(t2)
        assert "test_454" == list(t2)[0]

        shutil.rmtree("./tests/store")

    def test_describe(self):
        '''Function to test descriptors built into gblearn, held in gblearn/descriptors.py'''

        t1 = col("Test_1", "./tests/store")
        t1.read("./tests/test_data/ni.p455.out", 28,
                "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out',
                prefix="TEST")
        desc = "test"
        aid = "test_455"
        t1.describe(desc, rcut='huge', nmax='not_as_huge')
        res = t1.get(desc, aid, rcut='huge', nmax='not_as_huge')
        assert res == "test result"

        '''SOAP'''
        # FIXME make unit test that tests the functionality of SOAP
        t3 = col("Test_SOAP", "./tests/results")
        t3.read("./tests/test_data/ni.p455.out", 28,
                "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out')
        desc = 'soap'
        aid = '455'
        t3.describe(desc, rcut=5.0, nmax=9, lmax=9)
        fname = t3.store._generate_file_name(
            desc, aid, rcut=5.0, nmax=9, lmax=9)
        fpath = os.path.join(t3.store.root, desc, aid, fname)
        assert os.path.exists(fpath)
        #shutil.rmtree("./tests/results/")

        '''ASR'''
        # dummy SOAP array that ASR is run on np.array([[1,2,3,4],[3,4,5,6],[-1,0,4,2]])
        t4 = col("Test_ASR", "./tests/test_paths")
        t4.read("./tests/test_data/ni.p455.out", 28,
                "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out', prefix='fake')
        desc = 'asr'
        aid = 'fake_455'
        t4.describe(desc, needs_store=True, rcut=0, nmax=0, lmax=0)
        res = t4.get(desc, aid, rcut=0, nmax=0, lmax=0)
        exp = np.array([1, 2, 4, 4])  # expected  result of ASR
        assert np.array_equal(res, exp)
        fname = t4.store._generate_file_name(desc, aid, rcut=0, nmax=0, lmax=0)
        fpath = os.path.join(t4.store.root, desc, aid, fname)
        assert os.path.exists(fpath)
        shutil.rmtree("./tests/test_paths/asr")

        '''LER'''
        # FIXME make unit test that tests the functionality of LER
        t5 = col("Test_LER", "./tests/results")
        t5.read("./tests/test_data/ni.p455.out", 28, "lammps-dump-text",
            rxid=r'ni.p(?P<gbid>\d+).out')
        desc = 'ler'
        aid1 = '455'
        t5.describe(desc, needs_store=True, collection=t5,
                    eps=0.025, rcut=5.0, nmax=9, lmax=9)
        fname = t5.store._generate_file_name(
            desc, aid1, collection=t5, eps=0.025, rcut=5.0, nmax=9, lmax=9)
        fpath = os.path.join(t5.store.root, desc, aid1, fname)
        assert os.path.exists(fpath)
        shutil.rmtree("./tests/results/")
        #shutil.rmtree("./tests/test_paths/ler")

    def test_get(self):
        my_col = col("A", "./tests/results")
        result = "Random test result"
        desc = "test"
        aid = "12"
        my_col.store.store(result, desc, aid, arg_a=1, arg_b=2)
        ret_val3 = my_col.get(desc, aid, arg_a=1, arg_b=2)
        assert ret_val3 == result
