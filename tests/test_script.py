from gblearn.collection import Collection
from gblearn.store import ResultStore


c = Collection("A", "../store")
c.read(["../homer/ni.p454.out", "../homer/ni.p453.out"], 28, "lammps-dump-text", 
    rxid=r'ni.p(?P<gbid>\d+).out')
c.describe("soap",  rcut=5.0, nmax=9, lmax=9)
print(c.get_descriptor("453", "soap", rcut=5.0, nmax=9,lmax=9))


