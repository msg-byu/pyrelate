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
        t1 = col("t1", "./tests/store")
        rxid = r'ni.p(?P<gbid>\d+).out'
        import re
        c_rxid = re.compile(rxid)
        fn1 = "../homer/ni.p454.out"
        a1 = t1._read_aid(fn1, c_rxid)
        assert a1 == "454"
        shutil.rmtree("./tests/store")

    def test_read_aid_with_prefix(self):
        t1 = col("t1", "./tests/store")
        rxid = r'ni.p(?P<gbid>\d+).out'
        import re
        c_rxid = re.compile(rxid)
        fn1 = "../homer/ni.p454.out"
        prefix = "Pre"
        a2 = t1._read_aid(fn1, c_rxid, prefix)
        assert a2 == "pre_454"
        shutil.rmtree("./tests/store")

    def test_read_aid_no_regex(self):
        t1 = col("t1", "./tests/store")
        fn1 = "../homer/ni.p454.out"
        a4 = t1._read_aid(fn1, None)
        assert a4 == "ni.p454.out"
        shutil.rmtree("./tests/store")

    def test_read_aid_no_regex_with_prefix(self):
        t1 = col("t1", "./tests/store")
        fn1 = "../homer/ni.p454.out"
        prefix = "Test"
        a5 = t1._read_aid(fn1, None, prefix)
        assert a5 == "test_ni.p454.out"
        shutil.rmtree("./tests/store")

    def test_read_aid_invalid_regex(self):
        t1 = col("t1", "./tests/store")
        fn1 = "../homer/ni.p454.out"
        prefix = "Test"
        output2 = io.StringIO()
        sys.stdout = output2
        import re
        c_rxid3 = re.compile(r'ni.p(P<gbid>\d+).out')
        a6 = t1._read_aid(fn1, c_rxid3, prefix)
        assert "Regex found no pattern. Resolving to filename as aid.\n" == output2.getvalue()
        assert a6 == "test_ni.p454.out"
        shutil.rmtree("./tests/store")

    def test_read_list(self):
        t1 = col("Test_1", "./tests/store")
        # list of input files
        t1.read(["./tests/test_data/ni.p454.out", "./tests/test_data/ni.p453.out"], 28,
                "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out',
                prefix="TEST")
        assert 2 == len(t1)
        assert "test_454" == list(t1)[0]
        shutil.rmtree("./tests/store")

    def test_read_list_with_atomic_num_list(self):
        t1 = col("Test_1", "./tests/store")
        # list of input files with atomic num list (corresponding to one file)
        t1.read(["./tests/test_data/ni.p454.out", "./tests/test_data/ni.p453.out"], [28,28],
                "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out',
                prefix="TEST")
        assert 2 == len(t1)
        assert "test_454" == list(t1)[0]
        shutil.rmtree("./tests/store")

    def test_read_single_file(self):
        t1 = col("Test_1", "./tests/store")
        # read single file
        t1.read("./tests/test_data/ni.p455.out", 28,
                "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out',
                prefix="TEST")
        assert 1 == len(t1)
        shutil.rmtree("./tests/store")

    def test_read_directory(self):
        t1 = col("Test_1", "./tests/store")
        # read directory
        t1.read("./tests/test_data/sub1/", 28,
                "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out',
                prefix="TEST")
        assert 2 == len(t1)
        shutil.rmtree("./tests/store")

    def test_read_empty_dir_with_file(self):
        t1 = col("Test_1", "./tests/store")
        # empty directory / directory and file
        t1.read(["./tests/test_data/ni.p456.out", "./tests/test_data/empty"], 28,
                "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out',
                prefix="TEST")
        assert 1 == len(t1)
        shutil.rmtree("./tests/store")

    def test_read_empty_list(self):
        t1 = col("Test_1", "./tests/store")
        # empty list
        t1.read([], 28,
                "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out',
                prefix="TEST")
        assert 0 == len(t1)
        shutil.rmtree("./tests/store")

    def test_read_repeat_file(self):
        # will not read previously read file
        t1 = col("Test_1", "./tests/store")
        t1.read("./tests/test_data/ni.p456.out", 28,
                "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out',
                prefix="TEST")
        t1.read("./tests/test_data/ni.p456.out", 28,
                "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out',
                prefix="TEST")
        assert 1 == len(t1)
        shutil.rmtree("./tests/store")

    def test_read_nonex_directory(self):
        # non-existent directory
        t1 = col("Test_1", "./tests/store")
        output = io.StringIO()
        sys.stdout = output
        t1.read("definitely_wrong", 28,
                "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out',
                prefix="TEST")
        assert "Invalid file path, definitely_wrong was not read.\n" == output.getvalue()
        shutil.rmtree("./tests/store")

    def test_read_ASE_read_error(self):
        t1 = col("Test_1", "./tests/store")
        # ASE io.read() cannot automatically determine filetype
        self.assertRaises(StopIteration, col.read, t1, "./tests/test_data/ni.p457.out",
                          28, rxid=r'ni.p(?P<gbid>\d+).out', prefix="TEST")
        # missing parameter
        self.assertRaises(StopIteration, col.read, t1, "./tests/test_data/ni.p456.out",
                          "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out', prefix="TEST")
        shutil.rmtree("./tests/store")

    def test_read_no_filename(self):
        pass
        #test with xyz file

    def test_describe(self):
        t1 = col("Test_1", "./tests/store")
        t1.read("./tests/test_data/ni.p455.out", 28,
                "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out',
                prefix="TEST")
        desc = "test"
        aid = "test_455"
        t1.describe(desc, rcut='huge', nmax='not_as_huge')
        res = t1.get(desc, aid, rcut='huge', nmax='not_as_huge')
        assert res == "test result"
        shutil.rmtree("./tests/store")

    def test_describe_own_function(self):
        t1 = col("Test_1", "./tests/store")
        t1.read("./tests/test_data/ni.p455.out", 28,
                "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out', prefix="TEST")
        desc = "fake"
        aid = "test_455"
        from own_descriptor import Desc
        t1.describe(desc, Desc.fake_descriptor, arg1=1, arg2=2, arg3=3)
        res = t1.get(desc,aid, arg1=1, arg2=2, arg3=3)
        assert res == [1,2,3]
        shutil.rmtree("./tests/store")

    def test_get(self):
        my_col = col("A", "./tests/results")
        result = "Random test result"
        desc = "test"
        aid = "12"
        my_col.store.store(result, desc, aid, arg_a=1, arg_b=2)
        ret_val3 = my_col.get(desc, aid, arg_a=1, arg_b=2)
        assert ret_val3 == result
        shutil.rmtree("./tests/results")
