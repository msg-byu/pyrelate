import pyrelate.elements as elements
import numpy as np

class TestDescriptors():
    def test_seed_creation(self):
        elem = "Ni"
        a = elements.atoms(elem)
        assert np.all(a.get_pbc()) == True

    def test_no_lattice_data(self):
        elem = "Fe"
        try:
            a = elements.atoms(elem)
            assert False, f"Exception should be raised because {elem} is not hard-coded into elements.py"
        except:
            assert True




