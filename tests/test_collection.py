"""Tests collection.py
"""
from pyrelate.collection import AtomsCollection
import numpy as np
import sys
import shutil
import io
import re
import unittest
import os

'''Functions to help in writing and designing clear, functional unit tests'''


def _delete_store(my_col):
    '''Function to delete store generated in testing'''
    shutil.rmtree(my_col.store.root)


def _initialize_collection_and_read(aids):
    '''Initialize collection and read specified atoms files

    Parameters:
        aids (list of str): list of aid's for all ASE Atoms objects to be read into collection from  test_data
    '''
    my_col = AtomsCollection("Test", store="tests/store")
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
            assert my_col.get_description(aid, d, **kwargs) is not None
    return my_col


def _swap_x_y(positions):
    """Function to swap x and y coordinates of position"""
    transposed = positions.T.copy()
    temp = transposed[0].copy()
    transposed[0] = transposed[1].copy()
    transposed[1] = temp

    swapped_positions = transposed.T
    return swapped_positions


'''Toy descriptor functions to help in funtionality testing'''


def _test_descriptor(atoms, num=0, **kwargs):
    if num == 0:
        return 'test result 1', {}
    else:
        return 'test result 2', {}


def _processing_method(collection, based_on, method_name, **kwargs):
    # process collection of results
    new_string = method_name + "__"
    for aid in collection.aids():
        res = collection.get_description(aid, based_on[0], **based_on[1])  # str 'test result 1'
        new_string += res
        new_string += "_"
    info = {}
    return new_string, info


'''Unit tests'''


class TestCollection(unittest.TestCase):
    def test_subset_defaults(self):
        aids = ['454', '455']
        my_col = _initialize_collection_and_read(aids)
        data_loc = "tests/test_data/sub1"
        my_col.read(data_loc, 28, 'lammps-dump-text', rxid=r'ni.p(?P<aid>\d+).out')
        aids.append('456')

        new_col = my_col.subset(aids[:2])
        assert new_col.aids() == aids[:2]
        assert new_col.store.root == my_col.store.root
        assert new_col.name == my_col.name
        _delete_store(my_col)

    def test_subset_new_name(self):
        aids = ['454', '455']
        my_col = _initialize_collection_and_read(aids)

        new_col = my_col.subset(aids[:1], name="Test_sub")
        assert new_col.name == "test_sub"
        assert my_col.name == "test"
        assert new_col.aids() == aids[:1]
        _delete_store(my_col)

    def test_subset_new_store(self):
        aids = ['454', '455']
        my_col = _initialize_collection_and_read(aids)

        new_col = my_col.subset(aids[:1], store="tests/store_2")
        assert my_col.store.root != new_col.store.root
        assert new_col.aids() == aids[:1]
        _delete_store(my_col)
        _delete_store(new_col)

    def test_read_aid(self):
        '''Test _read_aid function'''
        my_col = AtomsCollection("Test", store="./tests/store")
        rxid = r'ni.p(?P<gbid>\d+).out'
        c_rxid = re.compile(rxid)
        filename = "./tests/test_data/ni.p454.out"
        aid = my_col._read_aid(filename, c_rxid)
        assert aid == "454"
        _delete_store(my_col)

    def test_read_aid_with_prefix(self):
        '''Test _read_aid, with prefix'''
        my_col = AtomsCollection("Test", store="./tests/store")
        rxid = r'ni.p(?P<gbid>\d+).out'
        c_rxid = re.compile(rxid)
        filename = "./tests/test_data/ni.p454.out"
        prefix = "Pre"
        aid = my_col._read_aid(filename, c_rxid, prefix)
        assert aid == "pre_454"
        _delete_store(my_col)

    def test_read_aid_no_regex(self):
        '''Test _read_aid, no regex'''
        my_col = AtomsCollection("Test", store="./tests/store")
        filename = "./tests/test_data/ni.p454.out"
        aid = my_col._read_aid(filename, None)
        assert aid == "ni.p454.out"
        _delete_store(my_col)

    def test_read_aid_no_regex_with_prefix(self):
        '''Test _read_aid, no regex but with prefix'''
        my_col = AtomsCollection("Test", store="./tests/store")
        filename = "./tests/test_data/ni.p454.out"
        prefix = "Test"
        aid = my_col._read_aid(filename, None, prefix)
        assert aid == "test_ni.p454.out"
        _delete_store(my_col)

    def test_read_aid_invalid_regex(self):
        '''Test _read_aid, invaid regex prints error and sets aid as filename'''
        my_col = AtomsCollection("Test", store="./tests/store")
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
        my_col = AtomsCollection("Test", store="./tests/store")
        my_col.read(["./tests/test_data/ni.p454.out", "./tests/test_data/ni.p455.out"], 28, rxid=r'ni.p(?P<gbid>\d+).out', prefix="TEST")
        assert 2 == len(my_col)
        assert "test_454" == list(my_col)[0]
        _delete_store(my_col)

    def test_read_list_with_atomic_num_list(self):
        '''Test read, list of input files with atomic number list'''
        my_col = AtomsCollection("Test", store="./tests/store")
        my_col.read(["./tests/test_data/ni.p454.out", "./tests/test_data/ni.p455.out"], [28, 28], rxid=r'ni.p(?P<gbid>\d+).out', prefix="TEST")
        assert 2 == len(my_col)
        assert "test_454" == list(my_col)[0]
        _delete_store(my_col)

    def test_read_single_file(self):
        '''Test read function, read single file'''
        my_col = AtomsCollection("Test", store="./tests/store")
        my_col.read("./tests/test_data/ni.p455.out", 28, rxid=r'ni.p(?P<gbid>\d+).out', prefix="TEST")
        assert 1 == len(my_col)
        _delete_store(my_col)

    def test_read_directory(self):
        '''Test read, read all input files in directory'''
        my_col = AtomsCollection("Test", store="./tests/store")
        my_col.read("./tests/test_data/sub1/", 28, rxid=r'ni.p(?P<gbid>\d+).out', prefix="TEST")
        assert 1 == len(my_col)
        _delete_store(my_col)

    def test_read_empty_dir_with_file(self):
        '''Test read, read empty directory + single file'''
        my_col = AtomsCollection("Test", store="./tests/store")
        my_col.read(["./tests/test_data/ni.p455.out", "./tests/test_data/empty"], 28, rxid=r'ni.p(?P<gbid>\d+).out', prefix="TEST")
        assert 1 == len(my_col)
        _delete_store(my_col)

    def test_read_empty_list(self):
        '''Test read, empty list'''
        my_col = AtomsCollection("Test", store="./tests/store")
        my_col.read([], 28, rxid=r'ni.p(?P<gbid>\d+).out', prefix="TEST")
        assert 0 == len(my_col)
        _delete_store(my_col)

    def test_read_repeat_file(self):
        '''Test read, repeat file, will not read previously read file'''
        my_col = _initialize_collection_and_read(['454'])
        my_col.read("./tests/test_data/ni.p454.out", 28, rxid=r'ni.p(?P<gbid>\d+).out')
        assert 1 == len(my_col)
        _delete_store(my_col)

    def test_read_nonexistent_directory(self):
        '''Test read, try to read nonexistent directory and throw error'''
        my_col = AtomsCollection("Test", store="./tests/store")
        output = io.StringIO()
        sys.stdout = output
        my_col.read("definitely_wrong", 28, rxid=r'ni.p(?P<gbid>\d+).out', prefix="TEST")
        assert "Invalid file path, definitely_wrong was not read.\n" == output.getvalue()
        _delete_store(my_col)

    def test_trim_correct_trim(self):
        '''Test trim function, check that a) some atoms are trimmed, and b) no outside values kept in the Atoms object'''
        aid = '454'
        xdim = 0
        trim_val = 3
        pad_val = 3

        my_col = _initialize_collection_and_read([aid])
        pre_trim_size = len(my_col[aid])
        my_col.trim(trim=trim_val, dim=xdim, pad=pad_val)
        post_trim_size = len(my_col[aid])
        assert pre_trim_size > post_trim_size, "No atoms were trimmed"

        positions = my_col[aid].get_positions()[:, xdim]
        my_col[aid].get_array("mask")
        for idx, atom in enumerate(positions):
            if(positions[idx] > (trim_val + pad_val) or positions[idx] < (-1 * (trim_val + pad_val))):
                assert False, "Atoms object not trimmed correctly"
        _delete_store(my_col)

    def test_trim_correct_pad(self):
        """Test that any atoms with a '0' value in the mask are supposed to be in the pad"""
        aid = '454'
        xdim = 0
        trim_val = 3
        pad_val = 3

        my_col = _initialize_collection_and_read([aid])
        my_col.trim(trim=trim_val, dim=xdim, pad=pad_val)

        positions = my_col[aid].get_positions()[:, xdim]
        mask = my_col[aid].get_array("mask")
        for idx, atom in enumerate(positions):
            if(mask[idx] == 0):
                if(positions[idx] < trim_val and positions[idx] > (trim_val * -1)):
                    assert False, "Mask was applied to atoms supposed to be included in final values"
        _delete_store(my_col)

    def test_trim_pad_True(self):
        '''Test that when pad=True, default is the expected value (equal to trim)'''
        aid = '454'
        xdim = 0
        trim_val = 3
        expected_pad_val = 3

        my_col = _initialize_collection_and_read([aid])
        len(my_col[aid])
        my_col.trim(trim=trim_val, dim=xdim, pad=True)

        positions = my_col[aid].get_positions()[:, xdim]
        # I dont know why this was here "mask = my_col[aid].get_array("mask")"
        for idx, atom in enumerate(positions):
            if(positions[idx] > (trim_val + expected_pad_val) or positions[idx] < (-1 * (trim_val + expected_pad_val))):
                assert False, "Pad not set to same as trim value"
        _delete_store(my_col)

    def test_trim_pad_False(self):
        '''Test that when pad=False, there is no pad'''
        aid = '454'
        xdim = 0
        trim_val = 3

        my_col = _initialize_collection_and_read([aid])
        my_col.trim(trim=trim_val, dim=xdim, pad=False)

        mask = my_col[aid].get_array("mask")
        assert np.count_nonzero(mask) == len(my_col[aid]), "Padding atoms included in mask when not expected"
        _delete_store(my_col)

    def test_trim_fail_invalid_trim(self):
        my_col = AtomsCollection("Test", "tests/store")
        # self.assertRaises(TypeError, AtomsCollection.trim, trim="string", dim=0)
        try:
            my_col.trim(trim="string", dim=0)
        except TypeError as e:
            assert e.__str__() == "Trim should be int or float type"
        else:
            assert False, "Expected type error not thrown"
        _delete_store(my_col)

    def test_trim_fail_invalid_pad(self):
        my_col = AtomsCollection("Test", store="tests/store")
        try:
            my_col.trim(trim=4, dim=0, pad="string")
        except TypeError as e:
            assert e.__str__() == "Pad should be int, float, or boolean type"
        else:
            assert False, "Expected type error not thrown"
        _delete_store(my_col)

    def test_trim_fail_invalid_dimension(self):
        aid = '454'
        my_col = _initialize_collection_and_read([aid])
        invalid_dim = 3
        try:
            my_col.trim(trim=4, dim=invalid_dim)
        except TypeError as e:
            assert e.__str__() == "Dimension should equal 0, 1, or 2"
        else:
            assert False, "Expected error not thrown"
        _delete_store(my_col)

    def test_trim_specify_diff_dimensions(self):
        """Test that specifying the dimension correctly trims different dimensions"""
        aid = '454'
        xdim = 0
        ydim = 1
        trim_val = 3

        my_col_A = _initialize_collection_and_read([aid])
        my_col_B = _initialize_collection_and_read([aid])
        new_positions = _swap_x_y(my_col_B[aid].get_positions())
        my_col_B[aid].set_positions(new_positions)

        my_col_A.trim(trim=trim_val, dim=xdim)
        my_col_B.trim(trim=trim_val, dim=ydim)

        mask_A = my_col_A[aid].get_array("mask")
        mask_B = my_col_B[aid].get_array("mask")
        assert np.array_equal(mask_A, mask_B), "Masks not equal, so atoms were trimmed differently for different dimensions"
        _delete_store(my_col_A)
        # deletes store b coincidentally

    def test_describe_own_function(self):
        '''Test using descriptor function not built into descriptors.py'''
        my_col = _initialize_collection_and_read(['455'])
        kwargs = {'num': 0, 'arg1': 1, 'arg2': 2, 'arg3': 3}
        my_col.describe('desc', fcn=_test_descriptor, **kwargs)
        res, info = my_col.get_description('455', 'desc', metadata=True, **kwargs)
        assert res == 'test result 1'
        assert info['desc_args'] == kwargs
        _delete_store(my_col)

    def test_describe_override(self):
        '''Put result in store, and check to make sure 'override' parameter overrides previous result'''
        kwargs = {'arg1': 1, 'arg2': 2, 'arg3': 3}
        my_col = _initialize_collection_and_read(['455'])
        my_col.store.store_description("fake result", {}, '455', "test", **kwargs)  # store result, can be overridden
        try:
            my_col.describe('test', fcn=_test_descriptor, override=True, **kwargs)
            res = my_col.get_description('455', 'test', **kwargs)
            assert res != "fake result"
            assert res == "test result 1"
        finally:
            _delete_store(my_col)

    def test_describe_trim_post_descriptor(self):
        aid = '455'
        my_col = _initialize_collection_and_read([aid])
        my_col.trim(trim=2, dim=0, pad=1)
        num_atoms_with_mask = len(my_col[aid])
        soapargs = {'rcut': 5.0, 'nmax': 3, 'lmax': 3}
        my_col.describe('soap', **soapargs)
        res, info = my_col.get_description(aid, 'soap', metadata=True, **soapargs)
        assert res is not None
        assert info['desc_args'] == soapargs
        assert num_atoms_with_mask > len(res), print("Result not correctly trimmed following description")
        _delete_store(my_col)

    def test_describe_specific_atomic_systems(self):
        a1 = '454'
        a2 = '455'
        my_col = _initialize_collection_and_read([a1, a2])
        kwargs = {'num': 0, 'arg1': 1, 'arg2': 2, 'arg3': 3}
        my_col.describe('desc', aid=a1, fcn=_test_descriptor, **kwargs)

        res1, info1 = my_col.get_description(a1, 'desc', metadata=True, **kwargs)
        assert res1 == 'test result 1'
        assert info1['desc_args'] == kwargs

        try:
            my_col.get_description(a2, 'desc', **kwargs)
            assert False, "No result was calculated previously so an exception should be raised"
        except FileNotFoundError:
            assert True

        _delete_store(my_col)

    def test_process(self):
        desc = "test"
        desc_args = {
            'num': 1
        }
        method = "my_method"
        method_args = {
            "a": 0,
            "b": 1,
            "method_name": method
        }
        my_col = _initialize_collection_and_describe([desc], ['454', '455'], **desc_args)
        res = my_col.process(method, (desc, desc_args), fcn=_processing_method, **method_args)
        assert res == "my_method__test result 2_test result 2_"

        _delete_store(my_col)

    def test_clear_method(self):
        '''Test clear, clear results created with method'''
        desc = "test"
        desc_args = {
            'num': 1
        }
        method = "my_method"
        method_args = {
            "a": 0,
            "b": 1,
            "method_name": method
        }
        my_col = _initialize_collection_and_describe([desc], ['454', '455'], **desc_args)
        my_col.process(method, (desc, desc_args), fcn=_processing_method, **method_args)
        my_col.clear(collection_name='Test', method=method)
        assert my_col.store.check_exists("Collections", 'Test', method) is False
        assert os.path.exists(os.path.join(my_col.store.root, "Collections", "Test", method)) is False
        _delete_store(my_col)

    def test_clear_description(self):
        '''Test clear, clear all results with specific description'''
        my_col = _initialize_collection_and_describe(['test'], ['454', '455'], arg1=1)
        my_col.clear(descriptor='test', aid='454')
        assert my_col.store.check_exists('Descriptions', '454', 'test') is False
        assert os.path.exists(os.path.join(my_col.store.root, 'Collections', '454', 'test')) is False
        _delete_store(my_col)

    def test_clear_specific_results_for_collection(self):
        '''Test clear, clear specific result from collection with specific collection_name, method, and based_on'''
        desc = "test"
        desc_args = {
            "num": 1
        }
        method = "my_method"
        method_args1 = {
            "a": 0,
            "b": 1,
            "method_name": method
        }
        method_args2 = {
            "a": 0,
            "b": 2,
            "method_name": method
        }
        my_col = _initialize_collection_and_describe([desc], ['454', '455'], **desc_args)
        my_col.process(method, (desc, desc_args), fcn=_processing_method, **method_args1)
        my_col.process(method, (desc, desc_args), fcn=_processing_method, **method_args2)
        my_col.clear(collection_name='test', method=method, based_on=(desc, desc_args), **method_args1)
        exists = my_col.store.check_exists("Collections", 'test', method, based_on=(desc, desc_args), **method_args1)
        assert exists is False
        exists = my_col.store.check_exists("Collections", 'test', method, based_on=(desc, desc_args), **method_args2)
        assert exists is True
        # check that can delete with default collection name
        my_col.clear(method=method, based_on=(desc, desc_args), **method_args2)
        exists = my_col.store.check_exists("Collections", 'test', method, based_on=(desc, desc_args), **method_args2)
        assert exists is False
        _delete_store(my_col)

    def test_clear_results_for_descriptor(self):
        '''Test clear, clear specific descriptor result with specific aid and descriptor'''
        my_col = _initialize_collection_and_describe(['test'], ['454'], arg1=1)
        my_col = _initialize_collection_and_describe(['test'], ['454'], arg2=2)
        my_col.clear(descriptor='test', aid='454', arg1=1)
        filename = my_col.store.check_exists('Descriptions', '454', 'test', arg1=1)
        assert filename is False
        filename = my_col.store.check_exists('Descriptions', '454', 'test', arg2=2)
        assert filename is True
        _delete_store(my_col)

    def test_clear_all(self):
        '''Test clear, clear all'''
        desc = "test"
        desc_args = {
            'num': 1
        }
        method = "my_method"
        method_args = {
            "a": 0,
            "b": 1,
            "method_name": method
        }
        my_col = _initialize_collection_and_describe([desc], ['454', '455'], **desc_args)
        my_col.process(method, (desc, desc_args), fcn=_processing_method, **method_args)
        my_col.clear_all()
        filename = my_col.store.check_exists("Collections", 'Test', method, based_on=desc, **method_args)
        assert filename is False
        filename = my_col.store.check_exists('Descriptions', '454', 'test', **desc_args)
        assert filename is False
        _delete_store(my_col)

    def test_aids(self):
        '''Test function to get list of all aid's in collection'''
        my_col = _initialize_collection_and_read(['454', '455'])
        aid_list = my_col.aids()
        assert type(aid_list) == list
        assert aid_list[0] == "454"
        assert len(aid_list) == 2
        _delete_store(my_col)
