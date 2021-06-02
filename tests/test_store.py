from pyrelate.store import Store
from pyrelate.collection import AtomsCollection
import os
import shutil
import unittest
# TODO finish updating unit tests


def _delete_store(store):
    shutil.rmtree(store.root)


def _test_descriptor(atoms, **kwargs):
    return 'test result'


def _initialize_collection_and_read(aids):
    '''Initialize collection and read specified atoms files'''
    my_col = AtomsCollection("Test", "tests/results")
    data_path = 'tests/test_data/ni.p{0:s}.out'
    for aid in aids:
        my_col.read(data_path.format(aid), 28, 'lammps-dump-text',
                    rxid=r'ni.p(?P<aid>\d+).out')
    return my_col


def _initialize_collection_and_describe(desc, aids, **kwargs):
    '''Initialize collection with specified descriptor and aids and describe'''
    my_col = _initialize_collection_and_read(aids)
    for d in desc:
        my_col.describe(d, fcn=_test_descriptor, **kwargs)
        for aid in aids:
            assert my_col.get_description(aid, d, **kwargs) is not None
    return my_col


class TestStore(unittest.TestCase):

    def test_init_store(self):
        '''Test initialization of store'''
        try:
            store = Store('./tests/results')
        except TypeError:
            assert False
        assert True
        _delete_store(store)

    def test_init_store_default(self):
        '''Test initializing store default'''
        try:
            store = Store()
        except TypeError:
            assert False
        assert True
        _delete_store(store)

    def test_init_store_expanduser(self):
        '''Test initializing store default'''
        try:
            store = Store("~/test_store")
        except TypeError:
            assert False
        assert True
        _delete_store(store)

    def test_based_on_is_correct_1(self):
        # if based_on is None, and not in dict, true
        store = Store("./tests/results")
        based_on = None
        info = {"other": 1}
        assert store._based_on_is_correct(based_on, info)
        _delete_store(store)

    def test_based_on_is_correct_2(self):
        # if based_on is None, and is in dict, false
        store = Store("./tests/results")
        based_on = None
        info = {"other": 1, "based_on_name": "name", "based_on_args": {"a": 1, "b": 2}}
        assert not store._based_on_is_correct(based_on, info)
        _delete_store(store)

    def test_based_on_is_correct_3(self):
        # if based_on is not None, and not in dict, false
        store = Store("./tests/results")
        desc_args = {"a": 1, "b": 2}
        based_on = ("desc", desc_args)
        info = {"other": 1}
        assert not store._based_on_is_correct(based_on, info)
        _delete_store(store)

    def test_based_on_is_correct_4(self):
        # if based_on is not None, and wrong based_on_args, false
        store = Store("./tests/results")
        desc_args = {"a": 1, "b": 2}
        name = "desc"
        based_on = (name, desc_args)
        info = {"other": 1, "based_on_name": name, "based_on_args": {"a": 1, "b": 22}}
        assert not store._based_on_is_correct(based_on, info)
        _delete_store(store)

    def test_based_on_is_correct_5(self):
        # if based_on is not None, and info matches, true
        store = Store("./tests/results")
        desc_args = {"a": 1, "b": 2}
        name = "desc"
        based_on = (name, desc_args)
        info = {"other": 1, "based_on_name": name, "based_on_args": desc_args}
        assert store._based_on_is_correct(based_on, info)
        _delete_store(store)

    def test_check_exists_description_true(self):
        store = Store("./tests/results")
        result = "Random test result"
        desc = "test_desc"
        aid = "111"
        kw1 = "option_1"
        kw2 = "option_2"
        info = {}
        store.store_description(result, info, aid, desc, a=kw1, b=kw2)

        assert store.check_exists("Descriptions", aid, desc, a=kw1, b=kw2)
        assert type(store.check_exists("Descriptions", aid, desc, a=kw1, b=kw2, explicit=True)) is str
        _delete_store(store)

    def test_check_exists_description_false(self):
        store = Store("./tests/results")
        desc = "test_desc"
        aid = "111"
        kw1 = "option_1"
        kw2 = "option_2"

        assert not store.check_exists("Descriptions", aid, desc, a=kw1, b=kw2)
        _delete_store(store)

    def test_check_exists_collection_true(self):
        store = Store("./tests/results")
        result = "Random test result"
        info = {
            "additional_info": [1, 2, 3, 4, 5]
        }
        desc = "test_desc"
        method = "test_method"
        name = "my_collection"
        desc_args = {
            "kw1": "option_1",
            "kw2": "option_2"
        }
        method_args = {
            "eps": 1,
            "num": 50
        }

        store.store_collection_result(result, info, method, name, (desc, desc_args), **method_args)

        assert store.check_exists("Collections", name, method, based_on=(desc, desc_args), **method_args)
        _delete_store(store)

    def test_check_exists_collection_false(self):
        store = Store("./tests/results")
        desc = "test_desc"
        method = "test_method"
        name = "my_collection"
        desc_args = {
            "kw1": "option_1",
            "kw2": "option_2"
        }
        method_args = {
            "eps": 1,
            "num": 50
        }

        assert not store.check_exists("Collections", name, method, based_on=(desc, desc_args), **method_args)
        _delete_store(store)

    def test_generate_default_file_name(self):
        '''Test _generate_default_file_name'''
        store = Store("./tests/results")
        aid = "455"
        desc = "soap"
        filename = store._generate_default_file_name(aid, desc)
        assert filename[:-26] == aid + "_" + desc
        assert filename[-4:] == ".pkl"
        _delete_store(store)

    def test_store_file(self):
        test = "thing"
        store = Store("./tests/results")
        path = os.path.join(store.root, "thing.pkl")
        store._store_file(test, path)
        assert os.path.exists(path)
        _delete_store(store)

    def test_store_description(self):
        # called in describe()
        # result, descriptor, aid, argmuments
        store = Store("./tests/results")
        result = "Random test result"
        desc = "test_desc"
        aid = "111"
        kw1 = "option_1"
        kw2 = "option_2"
        info = {}
        store.store_description(result, info, aid, desc, a=kw1, b=kw2)

        fpath = os.path.join(store.root, "Descriptions", aid, desc)
        os.path.join(store.root, "Descriptions", aid, desc)
        assert os.path.exists(fpath)

        directory = os.fsencode(fpath)
        for file in os.listdir(directory):
            filename = os.fsdecode(file)
            if filename[:-26] == aid + "_" + desc:
                fullpath = os.path.join(fpath, filename)
                assert os.path.exists(fullpath)
                break
        else:
            assert False, "No correct file found"
        _delete_store(store)

    def test_store_description_with_info(self):
        # called in describe()
        # result, descriptor, aid, argmuments
        store = Store("./tests/results")
        result = "Random test result"
        desc = "test_desc"
        aid = "111"
        kw1 = "option_1"
        kw2 = "option_2"
        info = {"num": 47, "important_info": 12, "fcn": _test_descriptor}
        store.store_description(result, info, aid, desc, a=kw1, b=kw2)

        fpath = os.path.join(store.root, "Descriptions", aid, desc)
        assert os.path.exists(fpath)

        directory = os.fsencode(fpath)
        for file in os.listdir(directory):
            filename = os.fsdecode(file)
            if filename[:-26] == aid + "_" + desc:
                fullpath = os.path.join(fpath, filename)
                assert os.path.exists(fullpath)
                info_fullpath = os.path.join(fpath, "info_" + filename)
                assert os.path.exists(info_fullpath)
                fetched_info = store._unpickle(info_fullpath)
                assert fetched_info['num'] == info['num']
                assert fetched_info['important_info'] == info['important_info']
                assert fetched_info['desc_args'] == {"a": kw1, "b": kw2}
                assert fetched_info['fcn'] == "_test_descriptor"
                break
        else:
            assert False, "No correct file found"
        _delete_store(store)

    def test_store_collection_description(self):
        # called in the process() method
        # result, info, collection name, arguments, descriptor_args
        store = Store("./tests/results")
        result = "Random test result"
        info = {
            "additional_info": [1, 2, 3, 4, 5]
        }
        desc = "test_desc"
        method = "test_method"
        name = "my_collection"
        desc_args = {
            "kw1": "option_1",
            "kw2": "option_2"
        }
        method_args = {
            "eps": 1,
            "num": 50
        }

        store.store_collection_result(result, info, method, name, (desc, desc_args), **method_args)

        fpath = os.path.join(store.root, "Collections", name, method)
        assert os.path.exists(fpath)

        directory = os.fsencode(fpath)
        for file in os.listdir(directory):
            filename = os.fsdecode(file)
            print(filename[:-26], name + "_" + method)
            if filename[:-26] == name + "_" + method:
                fullpath = os.path.join(fpath, filename)
                assert os.path.exists(fullpath)
                info_fullpath = os.path.join(fpath, "info_" + filename)
                assert os.path.exists(info_fullpath)
                fetched_info = store._unpickle(info_fullpath)
                assert fetched_info['method_args'] == method_args
                assert fetched_info['based_on_name'] == desc
                assert fetched_info['based_on_args'] == desc_args
                assert fetched_info['additional_info'] == info['additional_info']
                break
        else:
            assert False, "No correct file found"
        _delete_store(store)

    def test_store_collection_description_with_function_in_args(self):
        # called in the process() method
        # result, info, collection name, arguments, descriptor_args
        store = Store("./tests/results")
        result = "Random test result"
        info = {
            "additional_info": [1, 2, 3, 4, 5]
        }
        desc = "test_desc"
        method = "test_method"
        name = "my_collection"
        desc_args = {
            "kw1": "option_1",
            "kw2": "option_2"
        }
        method_args = {
            "eps": 1,
            "num": 50,
            "fcn": _test_descriptor
        }

        store.store_collection_result(result, info, method, name, (desc, desc_args), **method_args)

        fpath = os.path.join(store.root, "Collections", name, method)
        assert os.path.exists(fpath)

        directory = os.fsencode(fpath)
        for file in os.listdir(directory):
            filename = os.fsdecode(file)
            print(filename[:-26], name + "_" + method)
            if filename[:-26] == name + "_" + method:
                fullpath = os.path.join(fpath, filename)
                assert os.path.exists(fullpath)
                info_fullpath = os.path.join(fpath, "info_" + filename)
                assert os.path.exists(info_fullpath)
                fetched_info = store._unpickle(info_fullpath)
                assert fetched_info['method_args'] != method_args  # function converted to string of name
                assert fetched_info['method_args']['fcn'] == "_test_descriptor"
                assert fetched_info['based_on_name'] == desc
                assert fetched_info['based_on_args'] == desc_args
                assert fetched_info['additional_info'] == info['additional_info']
                break
        else:
            assert False, "No correct file found"
        _delete_store(store)

    def test_store_additional(self):
        pass

    def test_equal_args_true(self):
        dic1 = {"a": 1, "b": 2, "c": 3}
        dic2 = {"a": 1, "b": 2, "c": 3}
        store = Store("./tests/results")
        assert store._equal_args(dic1, dic2)
        _delete_store(store)

    def test_equal_args_false(self):
        dic1 = {"a": 1, "b": 2, "c": 3}
        dic2 = {"a": 1, "b": 2, "c": 4}
        store = Store("./tests/results")
        assert not store._equal_args(dic1, dic2)
        _delete_store(store)

    def test_equal_args_with_function(self):
        store = Store("./tests/results")
        func = _test_descriptor
        dic1 = {"func": func, "b": 2, "c": 3}
        dic2 = {"func": func.__name__, "b": 2, "c": 3}
        assert store._equal_args(dic1, dic2)
        _delete_store(store)

    def test_get_description(self):
        store = Store("./tests/results")
        result = "Random test result"
        desc = "test_desc"
        aid = "111"
        kwargs = {"a": "option_1", "b": "option_2"}
        info = {}
        store.store_description(result, info, aid, desc, **kwargs)

        res = store.get_description(aid, desc, **kwargs)

        assert res == result
        _delete_store(store)

    def test_get_description_metadata(self):
        store = Store("./tests/results")
        result = "Random test result"
        desc = "test_desc"
        aid = "111"
        kwargs = {"a": "option_1", "b": "option_2"}
        info = {}
        store.store_description(result, info, aid, desc, **kwargs)

        res, info = store.get_description(aid, desc, metadata=True, **kwargs)

        assert res == result
        assert info["desc_args"] == kwargs
        _delete_store(store)

    def test_get_collection_results(self):
        store = Store("./tests/results")
        result = "Random test result"
        info = {
            "additional_info": [1, 2, 3, 4, 5]
        }
        desc = "test_desc"
        method = "test_method"
        name = "my_collection"
        desc_args = {
            "kw1": "option_1",
            "kw2": "option_2"
        }
        method_args = {
            "eps": 1,
            "num": 50
        }

        store.store_collection_result(result, info, method, name, (desc, desc_args), **method_args)
        res = store.get_collection_result(method, name, (desc, desc_args), **method_args)

        assert res == result
        _delete_store(store)

    def test_get_collection_results_metadata(self):
        store = Store("./tests/results")
        result = "Random test result"
        info = {
            "additional_info": [1, 2, 3, 4, 5]
        }
        desc = "test_desc"
        method = "test_method"
        name = "my_collection"
        desc_args = {
            "kw1": "option_1",
            "kw2": "option_2"
        }
        method_args = {
            "eps": 1,
            "num": 50
        }

        store.store_collection_result(result, info, method, name, (desc, desc_args), **method_args)
        res, info = store.get_collection_result(method, name, (desc, desc_args), metadata=True, **method_args)

        assert res == result
        assert info["based_on_args"] == desc_args
        assert info["method_args"] == method_args
        _delete_store(store)

    def test_unpickle_path_and_fname(self):
        store = Store("./tests/results")
        test = "thing"
        fname = "thing.pkl"
        path = os.path.join(store.root, fname)
        store._store_file(test, path)

        fetched = store._unpickle(store.root, fname)
        assert fetched == test
        _delete_store(store)

    def test_unpickle_path(self):
        store = Store("./tests/results")
        test = "thing"
        fname = "thing.pkl"
        path = os.path.join(store.root, fname)
        store._store_file(test, path)

        fetched = store._unpickle(path)
        assert fetched == test
        _delete_store(store)

    def test_unpickle_does_not_exist(self):
        store = Store("./tests/results")
        try:
            store._unpickle("fakepath", "fake_fname")
        except FileNotFoundError:
            assert True
        else:
            assert False, "Expected error not thrown"
        finally:
            _delete_store(store)

    def test_unpickle_unpickling_error(self):
        import pickle
        store = Store("./tests/results/")
        fname = "fakepkl.pkl"
        try:
            store._unpickle("./tests", fname)
        except pickle.UnpicklingError:
            assert True
        else:
            assert False, "Expected error not thrown"
        finally:
            _delete_store(store)

    def test_clear_collection_result(self):
        '''Test clear, specific result'''
        store = Store("./tests/results")
        result1 = "Random test result1"
        result2 = "Random test result2"
        info = {
            "additional_info": [1, 2, 3, 4, 5]
        }
        desc = "test_desc"
        method = "test_method"
        name = "my_collection"
        desc_args = {
            "kw1": "option_1",
            "kw2": "option_2"
        }
        method_args1 = {
            "eps": 1,
            "num": 50
        }
        method_args2 = {
            "eps": 1,
            "num": 49
        }

        store.store_collection_result(result1, info, method, name, (desc, desc_args), **method_args1)
        store.store_collection_result(result2, info, method, name, (desc, desc_args), **method_args2)
        store.clear_collection_result(method, name, (desc, desc_args), **method_args1)

        try:
            store.get_collection_result(method, name, (desc, desc_args), **method_args1)
        except FileNotFoundError:
            assert True
        else:
            assert False, "Expected error not thrown"
        _delete_store(store)

    def test_clear_description_result(self):
        '''Test clear, clear all results for given discriptor and parameters for all aids in list'''
        store = Store("./tests/results")
        result = "Random test result"
        desc = "test_desc"
        aid = "111"
        kwargs1 = {"a": "option_1", "b": "option_2"}
        kwargs2 = {"a": "option_1", "b": "option_3"}
        info = {}
        store.store_description(result, info, aid, desc, **kwargs1)
        store.store_description(result, info, aid, desc, **kwargs2)
        store.clear_description_result(aid, desc, **kwargs1)

        try:
            store.get_description(aid, desc, **kwargs1)
        except FileNotFoundError:
            assert True
        else:
            assert False, "Expected error not thrown"
        _delete_store(store)

    def test_clear_method(self):
        '''Test clear, clear all results for given descriptor'''
        store = Store("./tests/results")
        result = "Random test result"
        info = {
            "additional_info": [1, 2, 3, 4, 5]
        }
        desc = "test_desc"
        method = "test_method"
        name = "my_collection"
        desc_args = {
            "kw1": "option_1",
            "kw2": "option_2"
        }
        method_args1 = {
            "eps": 1,
            "num": 50
        }
        method_args2 = {
            "eps": 1,
            "num": 49
        }
        store.store_collection_result(result, info, method, name, (desc, desc_args), **method_args1)
        store.store_collection_result(result, info, method, name, (desc, desc_args), **method_args2)

        store.clear_method(method, name)

        assert os.path.exists(os.path.join(store.root, "Collections", method, name)) is False
        _delete_store(store)

    def test_clear_description(self):
        '''Test clear, clear all results for given descriptor'''
        store = Store("./tests/results")
        result = "Random test result"
        desc = "test_desc"
        aid = "111"
        kwargs1 = {"a": "option_1", "b": "option_2"}
        kwargs2 = {"a": "option_1", "b": "option_3"}
        info = {}
        store.store_description(result, info, aid, desc, **kwargs1)
        store.store_description(result, info, aid, desc, **kwargs2)
        store.clear_description(aid, desc)

        assert os.path.exists(os.path.join(store.root, "Descriptions", aid, desc)) is False
        _delete_store(store)

    def test_clear_all(self):
        '''Test clear, clear all'''
        store = Store("./tests/results")
        result = "Random test result"
        info = {
            "additional_info": [1, 2, 3, 4, 5]
        }
        desc = "test_desc"
        method = "test_method"
        name = "my_collection"
        desc_args = {
            "kw1": "option_1",
            "kw2": "option_2"
        }
        method_args1 = {
            "eps": 1,
            "num": 50
        }

        aid = "111"
        kwargs1 = {"a": "option_1", "b": "option_2"}
        info = {}

        store.store_description(result, info, aid, desc, **kwargs1)
        store.store_collection_result(result, info, method, name, (desc, desc_args), **method_args1)
        store.clear_all()

        assert os.path.exists(os.path.join(store.root, "Descriptions", aid, desc)) is False
        assert os.path.exists(os.path.join(store.root, "Collections", method, name)) is False
        _delete_store(store)
