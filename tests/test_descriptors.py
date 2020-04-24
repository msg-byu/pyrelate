from gblearn.collection import AtomsCollection as col
import shutil
import os
import numpy as np

class TestDescriptors():
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
        shutil.rmtree("./tests/store")

    def test_soap(self):
        pass
        '''SOAP
        # FIXME make unit test that tests the functionality of SOAP
        t3 = col("Test_SOAP", "./tests/results")
        t3.read("./tests/test_data/ni.p455.out", 28,
                "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out')
        desc = 'soap'
        aid = '455'
        t3.describe(desc, rcut=5.0, nmax=9, lmax=9)
        fname = t3.store._generate_file_name(
            desc, aid, rcut=5.0, nmax=9, lmax=9)
        fpath = os.path.join(t3.store.root, desc, aid, fname)
        assert os.path.exists(fpath)
        shutil.rmtree("./tests/results/")
        '''

    def test_asr(self):
        '''ASR'''
        # dummy SOAP array that ASR is run on np.array([[1,2,3,4],[3,4,5,6],[-1,0,4,2]])
        t4 = col("Test_ASR", "./tests/test_paths")
        t4.read("./tests/test_data/ni.p455.out", 28,
                "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out', prefix='fake')
        desc = 'asr'
        aid = 'fake_455'
        t4.describe(desc, needs_store=True, rcut=0, nmax=0, lmax=0)
        res = t4.get(desc, aid, rcut=0, nmax=0, lmax=0)
        exp = np.array([1, 2, 4, 4])  # expected  result of ASR
        assert np.array_equal(res, exp)
        fname = t4.store._generate_file_name(desc, aid, rcut=0, nmax=0, lmax=0)
        fpath = os.path.join(t4.store.root, desc, aid, fname)
        assert os.path.exists(fpath)
        shutil.rmtree("./tests/test_paths/asr")

    def test_ler_runs(self):
        pass
        '''LER
        # FIXME make unit test that tests the functionality of LER
        t5 = col("Test_LER", "./tests/test_paths")
        t5.read("./tests/test_data/ni.p455.out", 28, "lammps-dump-text",
            rxid=r'ni.p(?P<gbid>\d+).out')
        desc = 'ler'
        aid1 = '455'
        t5.describe(desc, needs_store=True, collection=t5,
                    eps=0.025, rcut=5.0, nmax=9, lmax=9)
        fname = t5.store._generate_file_name(
            desc, aid1, collection=t5, eps=0.025, rcut=5.0, nmax=9, lmax=9)
        fpath = os.path.join(t5.store.root, desc, aid1, fname)
        assert os.path.exists(fpath)
        #shutil.rmtree("./tests/results/")
        shutil.rmtree("./tests/test_paths/ler")
        '''
    def test_ler_functionality(self):
        my_col = col("test_col", "./tests/results")
        my_col.read(["./tests/test_data/ni.p457.out","./tests/test_data/ni.p456.out"], 28, "lammps-dump-text", rxid=r'ni.p(?P<gbid>\d+).out', prefix="ler_unit_test")
        #store fake SOAP arrays
        rcut=0
        nmax=0
        lmax=0
        fake_mat1 = np.array([[1,2,1],[4,4,4],[5,4,5], [1,0,0]])
        fake_mat2 = np.array([[1,1,1],[10,10,9],[10,9,10],[1,1,2]])
        aid1="ler_unit_test_457"
        aid2="ler_unit_test_456"
        my_col.store.store(fake_mat1, "soap", aid1, rcut=rcut, nmax=nmax, lmax=lmax)
        my_col.store.store(fake_mat2, "soap", aid2, rcut=rcut, nmax=nmax, lmax=lmax)
        #compute LER
        seed =[0,0,0]
        my_col.describe("ler",needs_store=True, collection=my_col, eps=2.0, rcut=rcut,nmax=nmax,lmax=lmax, seed=seed)
        #check results
        U = my_col.get("ler", 'U', collection=my_col, eps=2.0, rcut=rcut, nmax=nmax, lmax=lmax, metric='euclidean', n_trees=10, search_k=-1)
        assert len(U['clusters']) == 4# 4 clusters
        print(len(U['centers']))
        ler1=my_col.get("ler", aid1,collection=my_col, eps=2.0, rcut=rcut, lmax=lmax, nmax=nmax, seed=seed)
        ler2=my_col.get("ler", aid2,collection=my_col, eps=2.0, rcut=rcut, lmax=lmax, nmax=nmax, seed=seed)
        #when sorting is implemented into LER these will be in a different order
        assert np.array_equal(ler1, np.array([1/4,1/2,1/4,0]))
        assert np.array_equal(ler2, np.array([1/2,0,0,1/2]))

        shutil.rmtree("./tests/results/")
