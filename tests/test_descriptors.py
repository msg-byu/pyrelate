from pyrelate.collection import AtomsCollection
import shutil
import os
import numpy as np


'''Functions to help in writing and designing clear, functional unit tests'''

def _delete_store(my_col):
    '''Function to delete store generated in testing'''
    shutil.rmtree(my_col.store.root)


def _initialize_collection_and_read(aids, store_loc="tests/store"):
    '''Initialize collection and read specified atoms files

    Parameters:
        aids (list of str): list of aid's for all ASE Atoms objects to be read into collection from  test_data
    '''
    my_col = AtomsCollection("Test", store_loc)
    data_path = 'tests/test_data/ni.p{0:s}.out'
    for aid in aids:
        my_col.read(data_path.format(aid), 28, 'lammps-dump-text', rxid=r'ni.p(?P<aid>\d+).out')
    my_col.trim(trim=2, dim=0, pad=1)
    return my_col

'''Unit Tests'''

class TestDescriptors():
    def test_soap(self):
        '''Test SOAP descriptor'''
        my_col = _initialize_collection_and_read(['455'])
        soapargs = {'rcut': 5.0, 'nmax': 9, 'lmax': 9}
        my_col.describe('soap', **soapargs)
        assert my_col.store.check_exists('soap', '455', **soapargs)
        res = my_col.get('soap', '455', **soapargs)
        assert type(res) is np.ndarray
        _delete_store(my_col)

    def test_asr(self):
        '''Test ASR descriptor'''
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
        '''Test ASR descriptor, norm_asr=True'''
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

    def test_sum(self):
        '''Test SUM descriptor'''
        my_col = _initialize_collection_and_read(['455'])
        soapargs = {'rcut': 0, 'nmax': 0, 'lmax': 0}
        fake_mat = np.array([[1, 2, 3, 4], [3, 4, 5, 6], [-1, 0, 4, 2]])
        my_col.store.store(fake_mat, "fake_soap", '455', **soapargs)
        my_col.describe('sum', res_needed='fake_soap', **soapargs)
        exp_res = np.array([3, 6, 12, 12])
        res = my_col.get('sum', '455', res_needed='fake_soap', **soapargs)
        assert np.array_equal(res, exp_res)
        _delete_store(my_col)

    def test_sum_fails(self):
        '''Test SUM descriptor'''
        my_col = _initialize_collection_and_read(['455'])
        soapargs = {'rcut': 0, 'nmax': 0, 'lmax': 0}
        my_col.describe('sum', res_needed='fake_soap', **soapargs)
        res = my_col.get('sum', '455', res_needed='fake_soap', **soapargs)
        assert res is None
        _delete_store(my_col)


    def test_ler_runs(self):
        '''Test LER function, previously computed SOAP'''
        #TODO store trimmed SOAP results to calculate LER on
        my_col = _initialize_collection_and_read(['455'], store_loc="tests/test_paths/") #has previously computed SOAP results stored here
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
        '''Test LER runs, check that using user specified SOAP function works'''
        my_col = _initialize_collection_and_read(
            ['455'], store_loc="tests/test_paths/")
        from pyrelate.descriptors import soap as soap_fcn
        lerargs = {
            'collection': my_col,
            'eps': 0.025,
            'rcut': 5.0,
            'nmax': 9,
            'lmax': 9,
            'soap_fcn': soap_fcn
        }
        my_col.describe('ler', **lerargs)
        assert my_col.store.check_exists('ler', '455', **lerargs)
        my_col.clear('ler')

    def test_ler_functionality(self):
        '''Test LER, see if gives expected results'''
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
