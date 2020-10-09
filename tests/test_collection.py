"""Tests collection.py
"""
from pyrelate.collection import AtomsCollection
import os.path
import numpy as np
import sys
import shutil
import io
import re
import unittest

'''Functions to help in writing and designing clear, functional unit tests'''


def _delete_store(my_col):
    '''Function to delete store generated in testing'''
    shutil.rmtree(my_col.store.root)


def _initialize_collection_and_read(aids):
    '''Initialize collection and read specified atoms files

    Parameters:
        aids (list of str): list of aid's for all ASE Atoms objects to be read into collection from  test_data
    '''
    my_col = AtomsCollection("Test", "tests/store")
    data_path = 'tests/test_data/ni.p{0:s}.out'
    for aid in aids:
        my_col.read(data_path.format(aid), 28, 'lammps-dump-text',
                    rxid=r'ni.p(?P<aid>\d+).out')
    return my_col


def _initialize_collection_and_describe(desc, aids, **kwargs):
    '''Initialize and describe collection with given aid's and descriptors

    Parameters:
        aids (list of str): list of aid's for all ASE Atoms objects to be read into collection from  test_data
        desc(list of str): list of all descriptors to be applied to collection
        kwargs(dict): "arguments" to be used in descriptor, only used in name of the results file
    '''
    my_col = _initialize_collection_and_read(aids)
    for d in desc:
        my_col.describe(d, fcn=_test_descriptor, **kwargs)
        for aid in aids:
            assert my_col.get(d, aid, **kwargs) != None
    return my_col


'''Toy descriptor functions to help in funtionality testing'''


def _test_descriptor(atoms, **kwargs):
    return 'test result'


def _test_descriptor_with_store(atoms, store, **kwargs):
    return 'another test result'


'''Unit tests'''


class TestCollection(unittest.TestCase):
    def test_read_aid(self):
        '''Test _read_aid function'''
        my_col = AtomsCollection("Test", "./tests/store")
        rxid = r'ni.p(?P<gbid>\d+).out'
        c_rxid = re.compile(rxid)
        filename = "./tests/test_data/ni.p454.out"
        aid = my_col._read_aid(filename, c_rxid)
        assert aid == "454"
        _delete_store(my_col)

    def test_read_aid_with_prefix(self):
        '''Test _read_aid, with prefix'''
        my_col = AtomsCollection("Test", "./tests/store")
        rxid = r'ni.p(?P<gbid>\d+).out'
        c_rxid = re.compile(rxid)
        filename = "./tests/test_data/ni.p454.out"
        prefix = "Pre"
        aid = my_col._read_aid(filename, c_rxid, prefix)
        assert aid == "pre_454"
        _delete_store(my_col)

    def test_read_aid_no_regex(self):
        '''Test _read_aid, no regex'''
        my_col = AtomsCollection("Test", "./tests/store")
        filename = "./tests/test_data/ni.p454.out"
        aid = my_col._read_aid(filename, None)
        assert aid == "ni.p454.out"
        _delete_store(my_col)

    def test_read_aid_no_regex_with_prefix(self):
        '''Test _read_aid, no regex but with prefix'''
        my_col = AtomsCollection("Test", "./tests/store")
        filename = "./tests/test_data/ni.p454.out"
        prefix = "Test"
        aid = my_col._read_aid(filename, None, prefix)
        assert aid == "test_ni.p454.out"
        _delete_store(my_col)

    def test_read_aid_invalid_regex(self):
        '''Test _read_aid, invaid regex prints error and sets aid as filename'''
        my_col = AtomsCollection("Test", "./tests/store")
        filename = "./tests/test_data/ni.p454.out"
        prefix = "Test"
        output = io.StringIO()
        sys.stdout = output
        invalid_rxid = r'ni.p(P<gbid>\d+).out'
        c_rxid = re.compile(invalid_rxid)
        aid = my_col._read_aid(filename, c_rxid, prefix)
        assert "Regex found no pattern. Resolving to filename as aid.\n" == output.getvalue()
        assert aid == "test_ni.p454.out"
        _delete_store(my_col)

    def test_read_list(self):
        '''Test read function, read list of input files'''
        my_col = AtomsCollection("Test", "./tests/store")
        my_col.read(["./tests/test_data/ni.p454.out", "./tests/test_data/ni.p455.out"], 28, "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out', prefix="TEST")
        assert 2 == len(my_col)
        assert "test_454" == list(my_col)[0]
        _delete_store(my_col)

    def test_read_list_with_atomic_num_list(self):
        '''Test read, list of input files with atomic number list'''
        my_col = AtomsCollection("Test", "./tests/store")
        my_col.read(["./tests/test_data/ni.p454.out", "./tests/test_data/ni.p455.out"], [28, 28], "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out', prefix="TEST")
        assert 2 == len(my_col)
        assert "test_454" == list(my_col)[0]
        _delete_store(my_col)

    def test_read_single_file(self):
        '''Test read function, read single file'''
        my_col = AtomsCollection("Test", "./tests/store")
        my_col.read("./tests/test_data/ni.p455.out", 28, "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out', prefix="TEST")
        assert 1 == len(my_col)
        _delete_store(my_col)

    def test_read_directory(self):
        '''Test read, read all input files in directory'''
        my_col = AtomsCollection("Test", "./tests/store")
        my_col.read("./tests/test_data/sub1/", 28, "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out', prefix="TEST")
        assert 1 == len(my_col)
        _delete_store(my_col)

    def test_read_empty_dir_with_file(self):
        '''Test read, read empty directory + single file'''
        my_col = AtomsCollection("Test", "./tests/store")
        my_col.read(["./tests/test_data/ni.p455.out", "./tests/test_data/empty"], 28, "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out', prefix="TEST")
        assert 1 == len(my_col)
        _delete_store(my_col)

    def test_read_empty_list(self):
        '''Test read, empty list'''
        my_col = AtomsCollection("Test", "./tests/store")
        my_col.read([], 28, "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out', prefix="TEST")
        assert 0 == len(my_col)
        _delete_store(my_col)

    def test_read_repeat_file(self):
        '''Test read, repeat file, will not read previously read file'''
        my_col = _initialize_collection_and_read(['454'])
        my_col.read("./tests/test_data/ni.p454.out", 28, "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out')
        assert 1 == len(my_col)
        _delete_store(my_col)

    def test_read_nonexistent_directory(self):
        '''Test read, try to read nonexistent directory and throw error'''
        my_col = AtomsCollection("Test", "./tests/store")
        output = io.StringIO()
        sys.stdout = output
        my_col.read("definitely_wrong", 28, "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out', prefix="TEST")
        assert "Invalid file path, definitely_wrong was not read.\n" == output.getvalue()
        _delete_store(my_col)

    def test_read_ASE_read_error(self):
        '''Test read, ASE read error if filetype not included '''
        my_col = AtomsCollection("Test", "./tests/store")
        # ASE io.read() cannot automatically determine filetype
        self.assertRaises(StopIteration, AtomsCollection.read, my_col, "./tests/test_data/ni.p454.out", 28, rxid=r'ni.p(?P<gbid>\d+).out', prefix="TEST")
        _delete_store(my_col)

    def test_read_no_filetype(self):
        pass  # TODO add unit test for automatic filetype reading through ASE
        # test with xyz file

    def test_describe_own_function(self):
        '''Test using descriptor function not built into descriptors.py'''
        my_col = _initialize_collection_and_read(['455'])
        kwargs = {'arg1': 1, 'arg2': 2, 'arg3': 3}
        my_col.describe('desc', fcn=_test_descriptor, **kwargs)
        res = my_col.get('desc', '455', **kwargs)
        assert res == 'test result'
        _delete_store(my_col)

    def test_describe_function_needs_store(self):
        '''Test that store is passed into descriptor function'''
        my_col = _initialize_collection_and_read(['455'])
        kwargs = {'arg1': 1, 'arg2': 2}
        try:
            my_col.describe('desc', fcn=_test_descriptor_with_store, **kwargs)
            #error thrown if not passed correctly
        finally:
            _delete_store(my_col)

    def test_describe_override(self):
        '''Put result in store, and check to make sure 'override' parameter overrides previous result'''
        kwargs = {'arg1': 1, 'arg2': 2, 'arg3': 3}
        my_col = _initialize_collection_and_describe( ['test'], ['455'], **kwargs)
        #assert my_col.get('test', '455', **kwargs) == "test result"
        # fcn name is not included in file name, so this will appear to be a previously computed result that can be overridden
        try:
            my_col.describe('test', fcn=_test_descriptor_with_store, override=True, **kwargs)
            assert my_col.get('test', '455', **kwargs) != "test result"
        finally:
            _delete_store(my_col)

    def test_clear_single_result(self):
        '''Test clear, clear single result with given descriptor, parameters, and aid'''
        my_col = _initialize_collection_and_describe(['test'], ['454', '455'], arg1=1)
        my_col.clear('test', '454', arg1=1)
        assert my_col.store.check_exists('test', '454', arg1=1) == False
        assert my_col.store.check_exists('test', '455', arg1=1) == True
        assert os.path.exists(os.path.join(my_col.store.root, 'test', '454')) == False
        assert os.path.exists(my_col.store.root) == True

    def test_clear_specific_results_for_collection(self):
        '''Test clear, no aid given, clears results of given parameters for all aids in collection'''
        my_col = _initialize_collection_and_describe(['test'], ['454', '455'], arg1=1)
        my_col.clear('test', arg1=1)
        assert my_col.store.check_exists('test', '454', arg1=1) == False
        assert my_col.store.check_exists('test', '455', arg1=1) == False
        assert os.path.exists(os.path.join( my_col.store.root, 'test', '454')) == False
        assert os.path.exists(os.path.join(my_col.store.root, 'test')) == False
        assert os.path.exists(my_col.store.root) == True

    def test_clear_results_for_descriptor(self):
        '''Test clear, clear all results for given descriptor'''
        my_col = _initialize_collection_and_describe(['test', 'test2'], ['454'], arg1=1)
        my_col.clear('test')
        assert my_col.store.check_exists('test', '454', arg1=1) == False
        assert my_col.store.check_exists('test2', '454', arg1=1) == True
        assert os.path.exists(os.path.join(my_col.store.root, 'test')) == False

    def test_clear_all(self):
        '''Test clear, clear all'''
        my_col = _initialize_collection_and_describe(['test', 'test2'], ['454', '455'], arg1=1)
        my_col.clear()
        assert os.path.exists(os.path.join(my_col.store.root, 'test', '454')) == False
        assert os.path.exists(os.path.join(my_col.store.root, 'test')) == False
        assert os.path.exists(os.path.join(my_col.store.root, 'test2')) == False
        assert os.path.exists(my_col.store.root) == True

    def test_get(self):
        '''Test AtomsCollection "get" function with aid provided'''
        my_col = _initialize_collection_and_describe(['test'], ['454'], a="arg1")
        try:
            fetched_result = my_col.get('test', '454', a="arg1")
            result = "test result"  # see function _test_descriptor
            assert fetched_result == result
        finally:
            _delete_store(my_col)

    def test_get_no_aid(self):
        '''Test AtomsCollection "get" function when no aid provided'''
        my_col = _initialize_collection_and_describe(['test', 'test2'], ['454', '455'], a="arg1")
        try:
            value = my_col.get('test', a="arg1")
            assert type(value) is dict
            assert len(value) is 2
        finally:
            _delete_store(my_col)

    def test_get_aid_not_str(self):
        my_col = _initialize_collection_and_describe(['test'], ['454'], a="arg1")
        try:
            self.assertRaises(ValueError, AtomsCollection.get, my_col, descriptor="test", idd=454, a="arg1")
        finally:
            _delete_store(my_col)

    def test_get_with_list(self):
        my_col = _initialize_collection_and_describe(['test'], ['454','455'], a="arg1")
        try:
            assert type(my_col.get('test', ['454','455'], a="arg1")) is dict
        finally:
            _delete_store(my_col)

    def test_aids(self):
        '''Test function to get list of all aid's in collection'''
        my_col = _initialize_collection_and_read(['454', '455'])
        aid_list = my_col.aids()
        assert type(aid_list) is list
        assert aid_list[0] == "454"
        assert len(aid_list) is 2
