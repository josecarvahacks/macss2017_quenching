# 02 May 2017 00:48:45
import h5py
import numpy as np
import matplotlib.pyplot as plt



def read_recdict_from_hdf5(h5file):
    """Read a dict of record arrays from hdf5."""
    f = h5py.File(h5file, "r")
    recdict = {}
    for grp, val in f.iteritems():
        print grp
        datasets = []
        dtypes = []
        for key in f[grp].keys():
            dset = f[grp][key][:]
            dtypename = f[grp][key].dtype.name
            dtype = (str(key), dtypename)
            datasets.append(dset)
            dtypes.append(dtype)
        print dtypes
        recdict[str(grp)] = np.rec.fromarrays(tuple(datasets), dtype=dtypes)
    f.close()
    return(recdict)

def read_mock(mockfile):
    """Read mock data into a numpy.recarray"""
    recdict = read_recdict_from_hdf5(mockfile)
    return(recdict['galaxy'])


def test_mock(mockfile):
    galrec = read_mock(mockfile)
    # get central galaxies
    iscen = galrec['lg_halo_mass'] > 1
    print np.sum(iscen)
    plt.hist(galrec['lg_halo_mass'][iscen], bins=100)
    plt.yscale('log')
    plt.show()
    plt.hist(galrec['lg_stellar_mass'][iscen], bins=100)
    plt.yscale('log')
    plt.show()

if __name__ == "__main__":
    mockfile = '/Users/ying/Dropbox/Public/iHODcatalog_bolshoi.h5'
    test_mock(mockfile)
