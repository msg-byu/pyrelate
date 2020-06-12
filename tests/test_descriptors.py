from pyrelate.collection import AtomsCollection
import shutil
import os
import numpy as np


def _delete_store(my_col):
    '''function to help in testing'''
    shutil.rmtree(my_col.store.root)


def _initialize_collection_and_read(aids, store_loc="tests/results"):
    '''initialize collection and read specified atoms files'''
    my_col = AtomsCollection("Test", store_loc)
    data_path = 'tests/test_data/ni.p{0:s}.out'
    for aid in aids:
        my_col.read(data_path.format(aid),
                    28, 'lammps-dump-text', rxid=r'ni.p(?P<aid>\d+).out')
    return my_col


class TestDescriptors():
    def test_soap(self):
        '''SOAP'''
        my_col = _initialize_collection_and_read(['455'])
        soapargs = {'rcut': 5.0, 'nmax': 9, 'lmax': 9}
        my_col.describe('soap', **soapargs)
        assert my_col.store.check_exists('soap', '455', **soapargs)
        res = my_col.get('soap', '455', **soapargs)
        assert type(res) is np.ndarray
        _delete_store(my_col)

    def test_asr(self):
        '''ASR'''
        my_col = _initialize_collection_and_read(['455'])
        soapargs = {'rcut': 0, 'nmax': 0, 'lmax': 0}
        fake_mat = np.array([[1, 2, 3, 4], [3, 4, 5, 6], [-1, 0, 4, 2]])
        my_col.store.store(fake_mat, "fake_soap", '455', **soapargs)
        my_col.describe('asr', res_needed='fake_soap', **soapargs)
        exp_res = np.array([1, 2, 4, 4])
        res = my_col.get('asr', '455', res_needed='fake_soap', **soapargs)
        assert np.array_equal(res, exp_res)
        _delete_store(my_col)

    def test_asr_normalize(self):
        '''ASR'''
        my_col = _initialize_collection_and_read(['455'])
        soapargs = {'rcut': 0, 'nmax': 0, 'lmax': 0}
        asrargs = {'res_needed': 'fake_soap', 'norm_asr': True}
        fake_mat = np.array([[1, 2, 3, 4], [3, 4, 5, 6], [-1, 0, 4, 2]])
        my_col.store.store(fake_mat, "fake_soap", '455', **soapargs)
        my_col.describe('asr', **asrargs, **soapargs)
        res = my_col.get('asr', '455', **asrargs, **soapargs)
        exp_mag = 6.0827625303  # np.sqrt(37) #Sqrt(4^2+4^2+2^2+1^1) = Sqrt(37)
        # [1,2,4,4] is the expected  result of ASR
        exp_res = np.array([1, 2, 4, 4]) / exp_mag
        assert np.all(np.isclose(res, exp_res))
        _delete_store(my_col)

    def test_ler_runs(self):
        '''LER'''
        my_col = _initialize_collection_and_read(
            ['455'], store_loc="tests/test_paths/")
        lerargs = {
            'collection': my_col,
            'eps': 0.025,
            'rcut': 5.0,
            'nmax': 9,
            'lmax': 9
        }
        my_col.describe('ler', **lerargs)
        assert my_col.store.check_exists('ler', '455', **lerargs)
        my_col.clear('ler')

    def test_ler_runs_pass_in_soapfcn(self):
        # LER
        my_col = _initialize_collection_and_read(
            ['455'], store_loc="tests/test_paths/")
        from pyrelate.descriptors import soap as soapfcn
        lerargs = {
            'collection': my_col,
            'eps': 0.025,
            'rcut': 5.0,
            'nmax': 9,
            'lmax': 9,
            'soapfcn': soapfcn
        }
        my_col.describe('ler', **lerargs)
        assert my_col.store.check_exists('ler', '455', **lerargs)
        my_col.clear('ler')

    def test_ler_functionality(self):
        my_col = _initialize_collection_and_read(['454', '455'])
        soapargs = {'rcut': 0, 'nmax': 0, 'lmax': 0}
        fake_mat1 = np.array([[1, 2, 1], [4, 4, 4], [5, 4, 5], [1, 0, 0]])
        fake_mat2 = np.array([[1, 1, 1], [10, 10, 9], [10, 9, 10], [1, 1, 2]])
        my_col.store.store(fake_mat1, "fake_soap", '454', **soapargs)
        my_col.store.store(fake_mat2, "fake_soap", '455', **soapargs)
        seed = [0, 0, 0]
        lerargs = {
            'collection': my_col,
            'eps': 2.0,
            'seed': seed,
            'res_needed': 'fake_soap'
        }
        my_col.describe("ler", **lerargs, **soapargs)
        ler1 = my_col.get("ler", '454', **lerargs, **soapargs)
        ler2 = my_col.get("ler", '455', **lerargs, **soapargs)
        assert len(ler1) == 4  # 4 clusters
        # when sorting is implemented into LER these will be in a different order
        assert np.array_equal(ler1, np.array([1 / 4, 1 / 2, 1 / 4, 0]))
        assert np.array_equal(ler2, np.array([1 / 2, 0, 0, 1 / 2]))
        _delete_store(my_col)
        # delete store?
