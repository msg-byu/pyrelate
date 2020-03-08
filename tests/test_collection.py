"""Tests collection.py
"""
from gblearn.collection import AtomsCollection as col
import os
import sys
import io
import unittest


class TestCollection(unittest.TestCase):
    def test_get_aid(self):
        '''
        Test get_aid function with different inputs (regex, prefix, filename),
            including incorrect Regex pattern
        '''
        t1 = col("t1")
        rxid = r'ni.p(?P<gbid>\d+).out'
        import re
        c_rxid = re.compile(rxid)
        fn1 = "../homer/ni.p454.out"
        a1 = t1._get_aid(fn1, c_rxid)
        assert a1 == "454"

        prefix = "Pre"
        a2 = t1._get_aid(fn1, c_rxid, prefix)
        assert a2 == "pre454"

        prefix = ""
        a3 = t1._get_aid(fn1, c_rxid, prefix)
        assert a3 == "454"

        c_rxid = None
        prefix = None
        a4 = t1._get_aid(fn1, c_rxid, prefix)
        assert a4 == "ni.p454.out"

        prefix = "Test_"
        a5 = t1._get_aid(fn1, c_rxid, prefix)
        assert a5 == "test_ni.p454.out"

        output2 = io.StringIO()
        sys.stdout = output2
        import re
        c_rxid2 = re.compile(r'ni.p(P<gbid>\d+).out')
        a6 = t1._get_aid(fn1, c_rxid2, prefix)
        assert "Regex found no pattern. Resolving to filename as aid.\n" == output2.getvalue()

    def test_read(self):
        '''
        Test read() function with various inputs types (list of file paths, directories, files, etc)
        '''
        t1 = col("Test_1")
        # list of input files
        t1.read(["./tests/test_data/ni.p454.out", "./tests/test_data/ni.p453.out"], 28,
                "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out',
                prefix="TEST")
        assert 2 == len(t1)
        assert "test454" == list(t1)[0]
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

    def test_describe(self):
        '''Function to test descriptors built into gblearn, held in gblearn/descriptors.py'''
        t1 = col("Test_1")
        t1.read("./tests/test_data/ni.p455.out", 28,
                "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out',
                prefix="TEST")
        res = t1.describe('test', rcut='huge', nmax='not_as_huge')
        assert "<class 'dict'>" == str(type(res))
        assert 'test result' == res['test455']
        assert 1 == len(res)

        '''SOAP'''
        t2 = col("Test_SOAP")
        t2.read("./tests/test_data/ni.p455.out", 28,
                "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out',
                prefix="")
        res = t1.describe('soap', rcut=5.0, nmax=9, lmax=9)
        assert "<class 'dict'>" == str(type(res))
