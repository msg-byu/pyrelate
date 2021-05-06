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
        # assert my_col.store.check_exists('soap', '455', **soapargs)
        res = my_col.get_description('455', 'soap', **soapargs)
        assert type(res) is np.ndarray
        _delete_store(my_col)

    def test_asr(self):
        '''Test ASR descriptor'''
        my_col = _initialize_collection_and_read(['455'])
        soapargs = {'rcut': 0, 'nmax': 0, 'lmax': 0}
        fake_mat = np.array([[1, 2, 3, 4], [3, 4, 5, 6], [-1, 0, 4, 2]])
        my_col.store.store_description(fake_mat, {}, '455', "fake_soap", **soapargs)
        my_col.process('asr', ('fake_soap', soapargs))
        exp_res = np.array([1, 2, 4, 4])
        res = my_col.get_collection_result('asr', ('fake_soap', soapargs))
        assert np.array_equal(res[0], exp_res)
        _delete_store(my_col)

    def test_asr_normalize(self):
        '''Test ASR descriptor, norm_asr=True'''
        my_col = _initialize_collection_and_read(['455'])
        soapargs = {'rcut': 0, 'nmax': 0, 'lmax': 0}
        asrargs = {'norm_asr': True}
        fake_mat = np.array([[1, 2, 3, 4], [3, 4, 5, 6], [-1, 0, 4, 2]])
        my_col.store.store_description(fake_mat, {}, '455', "fake_soap", **soapargs)
        my_col.process('asr', ('fake_soap', soapargs), **asrargs)
        res = my_col.get_collection_result('asr',('fake_soap', soapargs), **asrargs)
        exp_mag = 6.0827625303  # np.sqrt(37) #Sqrt(4^2+4^2+2^2+1^1) = Sqrt(37)
        # [1,2,4,4] is the expected  result of ASR
        exp_res = np.array([1, 2, 4, 4]) / exp_mag
        assert np.all(np.isclose(res[0], exp_res))
        _delete_store(my_col)

    def test_sum(self):
        '''Test SUM descriptor'''
        my_col = _initialize_collection_and_read(['455'])
        soapargs = {'rcut': 0, 'nmax': 0, 'lmax': 0}
        fake_mat = np.array([[1, 2, 3, 4], [3, 4, 5, 6], [-1, 0, 4, 2]])
        my_col.store.store_description(fake_mat, {}, "455", "fake_soap", **soapargs)
        my_col.process('sum', ('fake_soap', soapargs))
        exp_res = np.array([3, 6, 12, 12])
        res = my_col.get_collection_result('sum', ('fake_soap', soapargs))
        assert np.array_equal(res[0], exp_res)
        _delete_store(my_col)

    def test_sum_fails(self):
        '''Test SUM descriptor. When the result it is based on is not calculated, it should raise an error.'''
        my_col = _initialize_collection_and_read(['455'])
        soapargs = {'rcut': 0, 'nmax': 0, 'lmax': 0}
        try:
            my_col.process('sum', ('fake_soap', soapargs))
            assert False, "Exception should be raised."
        except FileNotFoundError:
            assert True
        _delete_store(my_col)


    def test_ler_functionality(self):
        '''Test LER, see if gives expected results'''
        my_col = _initialize_collection_and_read(['454', '455'])
        soapargs = {'rcut': 0, 'nmax': 0, 'lmax': 0}
        fake_mat1 = np.array([[-14, -13, -11], [4, 4, 4], [5, 4, 5], [1, 0, 1]])
        fake_mat2 = np.array([[1, 1, 1], [10, 10, 9], [10, 9, 10], [-14,-12,-12]])
        my_col.store.store_description(fake_mat1, {}, "454", "fake_soap", **soapargs)
        my_col.store.store_description(fake_mat2, {}, "455", "fake_soap", **soapargs)
        seed = [0, 0, 0]
        lerargs = {
            'eps': 0.3,
            'dissim_args':{"gamma":0.1},
            'seed': seed,
        }
        my_col.process("ler", ("fake_soap",soapargs), **lerargs)
        ler, info = my_col.get_collection_result("ler", ("fake_soap", soapargs), metadata=True, **lerargs)

        assert len(ler[0]) == 4  # 4 clusters
        assert info['num_clusters'] == 4
        # when sorting is implemented into LER these will be in a different order
        assert np.array_equal(ler[0], np.array([1/4, 1/2, 1/4, 0]))
        assert np.array_equal(ler[1], np.array([1/4, 0, 1/4, 1/2]))
        _delete_store(my_col)

    def test_ler_runs_pass_in_soapfcn(self):
        '''Test LER runs, check that using user specified SOAP function works'''
        my_col = _initialize_collection_and_read(['454', '455'])
        soapargs = {'rcut': 0, 'nmax': 0, 'lmax': 0}
        fake_mat1 = np.array([[-14, -13, -11], [4, 4, 4], [5, 4, 5], [1, 0, 1]])
        fake_mat2 = np.array([[1, 1, 1], [10, 10, 9], [10, 9, 10], [-14,-12,-12]])
        my_col.store.store_description(fake_mat1, {}, "454", "fake_soap", **soapargs)
        my_col.store.store_description(fake_mat2, {}, "455", "fake_soap", **soapargs)

        def soap_fcn(atoms, **kwargs):
            return [[0,0,0]]

        lerargs = {
            'eps': 0.3,
            'dissim_args': {"gamma" : 0.1},
            'soap_fcn':soap_fcn
        }
        try:
            my_col.process("ler", ("fake_soap",soapargs), **lerargs)
            ler, info = my_col.get_collection_result("ler", ("fake_soap", soapargs), metadata=True, **lerargs)
        finally:
            _delete_store(my_col)

        assert len(ler[0]) == 4  # 4 clusters
        assert info['num_clusters'] == 4
        # when sorting is implemented into LER these will be in a different order
        assert np.array_equal(ler[0], np.array([1/4, 1/2, 1/4, 0]))
        assert np.array_equal(ler[1], np.array([1/4, 0, 1/4, 1/2]))
