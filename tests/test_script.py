from collection import Collection
import descriptors as desc
import store
c = Collection("results", "../")
import inspect
#print(inspect.getfullargspec(desc.soap)[0])
#print(c._r_store.generate_file_name("soap", "555"))
c.read(["../homer/ni.p454.out", "../homer/ni.p453.out"], "lammps-dump-text",[28,28], 
    rxid=r'ni.p(?P<gbid>\d+).out', prefix="pre-")
c.describe("soap", "soap", rcut=5.0, nmax=9, lmax=9)

