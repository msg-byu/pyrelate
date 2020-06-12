"""Tests collection.py
"""
from pyrelate.collection import AtomsCollection
import os
import numpy as np
import sys
import shutil
import io
import unittest


def _delete_store(my_col):
    '''function to help in testing'''
    shutil.rmtree(my_col.store.root)


def _test_descriptor(atoms, **kwargs):
    return 'test result'


def _test_descriptor_with_store(atoms, store, **kwargs):
    return 'another test result'


def _initialize_collection_and_read(aids):
    '''initialize collection and read specified atoms files'''
    my_col = AtomsCollection("Test", "tests/results")
    data_path = 'tests/test_data/ni.p{0:s}.out'
    for aid in aids:
        my_col.read(data_path.format(aid),
                28, 'lammps-dump-text', rxid=r'ni.p(?P<aid>\d+).out')
    return my_col

def _initialize_collection_and_describe(desc, aids, **kwargs):
    '''initialize collection with specified descriptor and aids and describe with no args'''
    my_col = _initialize_collection_and_read(aids)
    my_col.describe(desc, fcn=_test_descriptor, **kwargs)
    for aid in aids:
        assert my_col.get(desc, aid, **kwargs) != None
    return my_col

class TestCollection(unittest.TestCase):
    '''
    def test_read_aid(self):
        t1 = AtomsCollection("t1", "./tests/store")
        rxid = r'ni.p(?P<gbid>\d+).out'
        import re
        c_rxid = re.compile(rxid)
        fn1 = "../homer/ni.p454.out"
        a1 = t1._read_aid(fn1, c_rxid)
        assert a1 == "454"
        shutil.rmtree("./tests/store")

    def test_read_aid_with_prefix(self):
        t1 = AtomsCollection("t1", "./tests/store")
        rxid = r'ni.p(?P<gbid>\d+).out'
        import re
        c_rxid = re.compile(rxid)
        fn1 = "../homer/ni.p454.out"
        prefix = "Pre"
        a2 = t1._read_aid(fn1, c_rxid, prefix)
        assert a2 == "pre_454"
        shutil.rmtree("./tests/store")

    def test_read_aid_no_regex(self):
        t1 = AtomsCollection("t1", "./tests/store")
        fn1 = "../homer/ni.p454.out"
        a4 = t1._read_aid(fn1, None)
        assert a4 == "ni.p454.out"
        shutil.rmtree("./tests/store")

    def test_read_aid_no_regex_with_prefix(self):
        t1 = AtomsCollection("t1", "./tests/store")
        fn1 = "../homer/ni.p454.out"
        prefix = "Test"
        a5 = t1._read_aid(fn1, None, prefix)
        assert a5 == "test_ni.p454.out"
        shutil.rmtree("./tests/store")

    def test_read_aid_invalid_regex(self):
        t1 = AtomsCollection("t1", "./tests/store")
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
        t1 = AtomsCollection("Test_1", "./tests/store")
        # list of input files
        t1.read(["./tests/test_data/ni.p454.out", "./tests/test_data/ni.p453.out"], 28,
                "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out',
                prefix="TEST")
        assert 2 == len(t1)
        assert "test_454" == list(t1)[0]
        shutil.rmtree("./tests/store")

    def test_read_list_with_atomic_num_list(self):
        t1 = AtomsCollection("Test_1", "./tests/store")
        # list of input files with atomic num list (corresponding to one file)
        t1.read(["./tests/test_data/ni.p454.out", "./tests/test_data/ni.p453.out"], [28,28],
                "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out',
                prefix="TEST")
        assert 2 == len(t1)
        assert "test_454" == list(t1)[0]
        shutil.rmtree("./tests/store")

    def test_read_single_file(self):
        t1 = AtomsCollection("Test_1", "./tests/store")
        # read single file
        t1.read("./tests/test_data/ni.p455.out", 28,
                "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out',
                prefix="TEST")
        assert 1 == len(t1)
        shutil.rmtree("./tests/store")

    def test_read_directory(self):
        t1 = AtomsCollection("Test_1", "./tests/store")
        # read directory
        t1.read("./tests/test_data/sub1/", 28,
                "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out',
                prefix="TEST")
        assert 2 == len(t1)
        shutil.rmtree("./tests/store")

    def test_read_empty_dir_with_file(self):
        t1 = AtomsCollection("Test_1", "./tests/store")
        # empty directory / directory and file
        t1.read(["./tests/test_data/ni.p456.out", "./tests/test_data/empty"], 28,
                "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out',
                prefix="TEST")
        assert 1 == len(t1)
        shutil.rmtree("./tests/store")

    def test_read_empty_list(self):
        t1 = AtomsCollection("Test_1", "./tests/store")
        # empty list
        t1.read([], 28,
                "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out',
                prefix="TEST")
        assert 0 == len(t1)
        shutil.rmtree("./tests/store")

    def test_read_repeat_file(self):
        # will not read previously read file
        t1 = AtomsCollection("Test_1", "./tests/store")
        t1.read("./tests/test_data/ni.p456.out", 28,
                "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out',
                prefix="TEST")
        t1.read("./tests/test_data/ni.p456.out", 28,
                "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out',
                prefix="TEST")
        assert 1 == len(t1)
        shutil.rmtree("./tests/store")

    def test_read_nonexistent_directory(self):
        # non-existent directory
        t1 = AtomsCollection("Test_1", "./tests/store")
        output = io.StringIO()
        sys.stdout = output
        t1.read("definitely_wrong", 28,
                "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out',
                prefix="TEST")
        assert "Invalid file path, definitely_wrong was not read.\n" == output.getvalue()
        shutil.rmtree("./tests/store")

    def test_read_ASE_read_error(self):
        t1 = AtomsCollection("Test_1", "./tests/store")
        # ASE io.read() cannot automatically determine filetype
        self.assertRaises(StopIteration, AtomsCollection.read, t1, "./tests/test_data/ni.p457.out",
                          28, rxid=r'ni.p(?P<gbid>\d+).out', prefix="TEST")
        # missing parameter
        self.assertRaises(StopIteration, AtomsCollection.read, t1, "./tests/test_data/ni.p456.out",
                          "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out', prefix="TEST")
        shutil.rmtree("./tests/store")

    def test_read_no_filetype(self):
        pass #TODO add unit test for automatic filetype reading through ASE
        #test with xyz file
    '''
    def test_describe_built_in_function(self):
        #FIXME change so test is not in descriptors.py, use one of functions already there
        #maybe get rid of??
        t1 = AtomsCollection("Test_1", "./tests/store")
        t1.read("./tests/test_data/ni.p455.out", 28,
                "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out',
                prefix="TEST")
        desc = "test"
        aid = "test_455"
        args = {
            'arg1' : 'huge',
            'arg2' : 'not_as_huge'
        }
        t1.describe(desc, **args)
        res = t1.get(desc, aid, **args)
        assert res == "test result"
        shutil.rmtree("./tests/store")

    def test_describe_own_function(self):
        '''Test using descriptor function not built into descriptors.py'''
        my_col = _initialize_collection_and_read(['455'])
        kwargs = {'arg1' : 1,'arg2' : 2,'arg3' : 3}
        my_col.describe('desc', fcn=_test_descriptor, **kwargs)
        res = my_col.get('desc', '455', **kwargs)
        assert res == 'test result'
        _delete_store(my_col)

    def test_describe_function_needs_store(self):
        '''Test that store is passed into descriptor function'''
        my_col = _initialize_collection_and_read(['455'])
        kwargs = {'arg1' : 1,'arg2' : 2}
        try:
            my_col.describe('desc', fcn=_test_descriptor_with_store, **kwargs)
            assert True
        except TypeError:#TypeError is thrown when 'store' parameter is not correctly included in descriptor function
            assert False
        _delete_store(my_col)

    def test_describe_override(self):
        '''Put result in store, and check to make sure override indeed overrides it'''
        kwargs = {'arg1' : 1,'arg2' : 2,'arg3' : 3}
        my_col = _initialize_collection_and_describe('test', ['455'], **kwargs)
        assert my_col.get('test', '455', **kwargs) == "test result"
        #fcn name is not included in file name, so this will appear to be a previously computed result that can be overridden
        my_col.describe('test', fcn=_test_descriptor_with_store, override=True, **kwargs)
        assert my_col.get('test', '455', **kwargs) != "test result"
        _delete_store(my_col)

    def test_clear_single_result(self):
        pass

    def test_clear_specific_results_for_collection(self):
        pass

    def test_clear_results_for_descriptor_type(self):
        pass

    def test_clear_all(self):
        pass

    def test_get(self):
        my_col = AtomsCollection("A", "./tests/results")
        result = "Random test result"
        desc = "test"
        aid = "12"
        my_col.store.store(result, desc, aid, arg_a=1, arg_b=2)
        ret_val3 = my_col.get(desc, aid, arg_a=1, arg_b=2)
        assert ret_val3 == result
        shutil.rmtree("./tests/results")

    def test_get_no_aid(self):
        my_col = AtomsCollection("A", "./tests/results")
        result = "Random test result"
        desc = "test"
        aid = "12"
        aid2 = "13"
        my_col["12"] = "test"
        my_col["13"] = "test2"
        my_col.store.store(result, desc, aid, arg_a=1, arg_b=2)
        my_col.store.store(result, desc, aid2, arg_a=1, arg_b=2)
        ret_val3 = my_col.get(desc, arg_a=1, arg_b=2)
        assert type(ret_val3) is dict
        assert len(ret_val3) is 2
        shutil.rmtree("./tests/results")

    def test_aids(self):
        t1 = AtomsCollection("Test_1", "./tests/test_paths")
        t1.read(["./tests/test_data/ni.p454.out", "./tests/test_data/ni.p453.out"], 28,
                "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out')
        aid_list = t1.aids()
        assert type(aid_list) is list
        assert int(aid_list[0]) == int("453")
        assert len(aid_list) is 2
