# 02 May 2017 00:48:45
import h5py
import numpy as np
import matplotlib.pyplot as plt
try:
    from zypy.zycosmo import get_halofuncs, get_cosmo, CosmoParams
    has_hmf = True
except ImportError:
    has_hmf = False
    print("Cosmology/HMF code is required.")

try:
    from hod.shmr import SHMR_Leauthaud, HSMR, HSMR_split
    from hod.predict import get_f_sigma_lnMs
    has_shmr = True
except ImportError:
    print("SHMR code is required.")
    has_shmr = False


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
    mockdata = recdict['galaxy']
    print "The columns are: "
    print mockdata.dtype.names
    return(mockdata)

def read_mock_hmf(mockfile, mmin=1.e9, mmax=1.e16, nmbin=101, h=0.701):
    galrec = read_mock(mockfile)
    iscen = galrec['lg_halo_mass'] > 1
    _Mh_arr = np.logspace(np.log10(mmin), np.log10(mmax), nmbin)
    Mh_arr = np.sqrt(_Mh_arr[1:] * _Mh_arr[:-1])
    dn_arr = np.histogram(galrec['lg_halo_mass'][iscen] - np.log10(h), bins=np.log10(_Mh_arr))[0]
    dndMh_arr = dn_arr / (_Mh_arr[1:] - _Mh_arr[:-1])
    return(Mh_arr, dndMh_arr)


def test_mock_hmf(mockfile):
    """Check the halo mass function in the mock."""
    galrec = read_mock(mockfile)
    # get central galaxies that correspond to main dark matter halos.
    iscen = galrec['lg_halo_mass'] > 1
    print 'total number of halos: %8d' % np.sum(iscen)
    # halo masses are in units of Msun/h
    plt.hist(galrec['lg_halo_mass'][iscen], bins=100, alpha=0.5)
    # plt.hist(galrec['lg_stellar_mass'][iscen], bins=100, alpha=0.5)
    if has_hmf:
        # compare with theory
        rcube = 250.0 # Mpc/h
        cosmo = CosmoParams(omega_M_0=0.27, sigma_8=0.82, h=0.70, omega_b_0=0.0469, n=0.95, set_flat=True)
        M_arr, dndM_arr = get_halofuncs(z=0.1, cosmo=cosmo, DELTA_HALO=200.0, mmin=1.e9, mmax=1.e16, nmbin=100)[:2]
        dlogm = np.log(M_arr[1]) - np.log(M_arr[0])
        nhalo = dndM_arr * M_arr * cosmo.h * dlogm * rcube**3 / (cosmo.h)**3
        plt.plot(np.log10(M_arr * cosmo.h), nhalo, 'r-')
    plt.yscale('log')
    plt.xlabel(r"$M_h\;[M_\odot/h]$")
    plt.ylabel(r"$N$")
    plt.ylim(1e0, 1e6)
    plt.show()


def test_mock_shmr(mockfile):
    """Check the stellar to halo mass relation in the mock."""
    galrec = read_mock(mockfile)
    iscen = galrec['lg_halo_mass'] > 1
    lgmh = galrec['lg_halo_mass'][iscen]
    lgms = galrec['lg_stellar_mass'][iscen]
    lgmh_bins = np.linspace(11.4, 15.0, 35)
    lgmh_cens = (lgmh_bins[1:] + lgmh_bins[:-1]) / 2.0
    lgms_cens = np.empty_like(lgmh_cens)
    lgms_scas = np.empty_like(lgmh_cens)
    for i in xrange(lgmh_cens.size):
        sel = (lgmh >= lgmh_bins[i]) & (lgmh < lgmh_bins[i+1])
        lgms_cens[i] = np.mean(lgms[sel])
        lgms_scas[i] = np.std(lgms[sel])
    if has_shmr:
        lgMs_0 = 10.30790
        lgMh_1 = 12.09899
        beta = 0.33272
        delta = 0.440
        gamma = 1.20579
        Mh_1 = 10**lgMh_1  # Msun rather than Msun/h
        Ms_0 = 10**lgMs_0  # this is always Msun/h^2
        shmr = SHMR_Leauthaud(Mh_1=Mh_1, Ms_0=Ms_0, beta=beta, delta=delta, gamma=gamma)
        # scatter in the SHMR
        sigma_lnMs = 0.49730
        eta = -0.04104
        lgMh_sca = lgMh_1
        Mh_sca = 10**lgMh_sca
        f_sigma_lnMs = get_f_sigma_lnMs(sigma_lnMs=sigma_lnMs, eta=eta, Mh_sca=Mh_sca)
        #
        mharr = 10**lgmh_cens  / 0.701
        lnmsarr = shmr.log_stellarmass_mean(np.log(mharr))
        lgmsarr = lnmsarr / np.log(10.0)
        lgmssca = f_sigma_lnMs(mharr) / np.log(10.0)
        plt.plot(lgmh_cens, lgmsarr, 'r-')
        plt.plot(lgmh_cens, lgmsarr + lgmssca , 'r--')
        plt.plot(lgmh_cens, lgmsarr - lgmssca , 'r--')
    plt.plot(lgmh_cens, lgms_cens, 'k-')
    plt.plot(lgmh_cens, lgms_cens+lgms_scas, 'k--')
    plt.plot(lgmh_cens, lgms_cens-lgms_scas, 'k--')
    plt.xlabel(r"$M_h\;[M_\odot/h]$")
    plt.ylabel(r"$M_*\;[M_\odot/h^2]$")
    plt.show()

def test_mock_hsmr(mockfile):
    """Check the halo to stellar mass relation in the mock."""
    galrec = read_mock(mockfile)
    iscen = galrec['lg_halo_mass'] > 1
    lgmh = galrec['lg_halo_mass'][iscen]
    lgms = galrec['lg_stellar_mass'][iscen]
    lgms_bins = np.linspace(9.8, 12.0, 31)
    lgms_cens = (lgms_bins[1:] + lgms_bins[:-1]) / 2.0
    lgmh_cens = np.empty_like(lgms_cens)
    lgmh_scas = np.empty_like(lgms_cens)
    for i in xrange(lgms_cens.size):
        sel = (lgms >= lgms_bins[i]) & (lgms < lgms_bins[i+1])
        lgmh_cens[i] = np.mean(lgmh[sel])
        lgmh_scas[i] = np.std(lgmh[sel])
    if has_shmr:
        lgMs_0 = 10.30790
        lgMh_1 = 12.09899
        beta = 0.33272
        delta = 0.440
        gamma = 1.20579
        Mh_1 = 10**lgMh_1  # Msun rather than Msun/h
        Ms_0 = 10**lgMs_0  # this is always Msun/h^2
        shmr = SHMR_Leauthaud(Mh_1=Mh_1, Ms_0=Ms_0, beta=beta, delta=delta, gamma=gamma)
        # scatter in the SHMR
        sigma_lnMs = 0.49730
        eta = -0.04104
        lgMh_sca = lgMh_1
        Mh_sca = 10**lgMh_sca
        f_sigma_lnMs = get_f_sigma_lnMs(sigma_lnMs=sigma_lnMs, eta=eta, Mh_sca=Mh_sca)
        #
        cosmo = CosmoParams(omega_M_0=0.27, sigma_8=0.82, h=0.70, omega_b_0=0.0469, n=0.95, set_flat=True)
        M_arr, dndM_arr = get_halofuncs(z=0.1, cosmo=cosmo, DELTA_HALO=200.0, mmin=1.e9, mmax=1.e16, nmbin=100)[:2]
        # HSMR
        hsmr = HSMR(shmr, f_sigma_lnMs, dndM_arr, M_arr, lgmsmin=8.0, lgmsmax=13.0, dlgms=0.02)
        lnMh_mean, lnMh_mean2, lnMh_med, sigma_lnMh_low, sigma_lnMh_upp = hsmr.get_plnMh_at_lnMs()
        #
        h = 0.701
        lgms_arr = hsmr.lnMs_arr / np.log(10.0)
        lgmharr = lnMh_mean / np.log(10.0) + np.log10(h)
        sigupp = sigma_lnMh_upp / np.log(10.0)
        siglow = sigma_lnMh_low / np.log(10.0)
        plt.plot(lgms_arr, lgmharr, 'r-')
        plt.plot(lgms_arr, lgmharr + sigupp , 'r--')
        plt.plot(lgms_arr, lgmharr - siglow , 'r--')
        #
    plt.plot(lgms_cens, lgmh_cens, 'k-')
    plt.plot(lgms_cens, lgmh_cens+lgmh_scas, 'k--')
    plt.plot(lgms_cens, lgmh_cens-lgmh_scas, 'k--')
    plt.xlabel(r"$M_*\;[M_\odot/h^2]$")
    plt.ylabel(r"$M_h\;[M_\odot/h]$")
    plt.show()


if __name__ == "__main__":
    mockfile = '/Users/ying/Dropbox/Public/iHODcatalog_bolshoi.h5'
    # test_mock_hmf(mockfile)
    test_mock_shmr(mockfile)
    # test_mock_hsmr(mockfile)
