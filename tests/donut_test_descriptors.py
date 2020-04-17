from gblearn.collection import AtomsCollection as col
import shutil
import os
import numpy as np

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

        '''SOAP'''
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
        #shutil.rmtree("./tests/results/")

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

        '''LER'''
        # FIXME make unit test that tests the functionality of LER
        t5 = col("Test_LER", "./tests/results")
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
        shutil.rmtree("./tests/results/")
        #shutil.rmtree("./tests/test_paths/ler")
