"""Tests collection.py
"""
from gblearn.collection import AtomsCollection as col
import os
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
        Test read() function with various inputs types (list of file paths, directories, files, etc)
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
        #TODO test for a list of atomic numbers
        t2 = col("Test_2", "./tests/store")
        t2.read(["./tests/test_data/ni.p454.out", "./tests/test_data/ni.p453.out"], [28,28],
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
        '''
        #can I have a cooked up example for this one?
        #test it does something
        t2 = col("Test_SOAP", "./tests/store")
        t2.read("./tests/test_data/ni.p455.out", 28,
                "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out',
                prefix="")
        res = t1.describe('soap', rcut=5.0, nmax=9, lmax=9)
        assert "<class 'dict'>" == str(type(res))

	    #test running full descriptor and saving to file (and deleting after)
        t3 = col("Test_SOAP_2")
        rs = ResultStore("./tests/results/")
        t3.read("./tests/test_data/ni.p455.out", 28,
                "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out')
        desc = 'soap'
        aid = '455'
        res = t3.describe(desc ,rs, rcut=5.0, nmax=9, lmax=9)
        fname = rs.generate_file_name(desc, aid, rcut=5.0, nmax=9, lmax=9)
        fname = fname + '.npy'
        fpath = os.path.join(rs.root, desc, aid, fname)
        assert os.path.exists(fpath)
        shutil.rmtree("./tests/results/")
        '''
        '''ASR'''
        #cook up random soap vector, and compute asr by hand


    def test_get(self):
        my_col = col("A", "./tests/results")
        result = "Random test result"
        desc = "test"
        aid = "12"
        my_col.store.store(result, desc, aid, arg_a=1, arg_b=2)
        ret_val3 = my_col.get(desc, aid, arg_a=1, arg_b=2)
        assert ret_val3 == result
