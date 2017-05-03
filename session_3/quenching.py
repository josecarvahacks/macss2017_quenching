# 02 May 2017 16:39:54
import sys
sys.path.insert(1, '../session_2/')
try:
    from hod.fred import RedFractionGenerator
    from zypy.zycosmo import get_halofuncs, get_cosmo, CosmoParams
    from hod.shmr import SHMR_Leauthaud, HSMR, HSMR_split
    from hod.predict import get_f_sigma_lnMs
    has_fred = True
except ImportError:
    has_fred = False
import numpy as np
import matplotlib.pyplot as plt
from read_mock import read_mock, read_mock_hmf



"""Implement Halo Quenching."""


def is_red(lgms, gcolor):
    cut = 0.8*(lgms/10.5)**0.6
    isred = gcolor >= cut
    return(isred)

def test_mock_hsmr_split(mockfile):
    """Check the color-split version of halo to stellar mass relations in the mock."""
    galrec = read_mock(mockfile)
    iscen = galrec['lg_halo_mass'] > 0
    lgmh = galrec['lg_halo_mass'][iscen]
    lgms = galrec['lg_stellar_mass'][iscen]
    gcolor = galrec['g-r'][iscen]
    lgms_bins = np.linspace(9.8, 12.5, 20)
    # lgms_bins = np.array([9.5, 9.8, 10.0, 10.2, 10.4, 10.6, 10.8, 11.0, 11.2, 11.4, 12.0])
    lgms_cens = (lgms_bins[1:] + lgms_bins[:-1]) / 2.0
    # initiate arrays
    lgms_cens_red = np.empty_like(lgms_cens)
    lgmh_cens_red = np.empty_like(lgms_cens)
    lgmh_err_red = np.empty_like(lgms_cens)
    lgms_cens_blue = np.empty_like(lgms_cens)
    lgmh_cens_blue = np.empty_like(lgms_cens)
    lgmh_err_blue = np.empty_like(lgms_cens)
    # go thru each stellar mass bin
    for i in xrange(lgms_cens.size):
        sel = (lgms >= lgms_bins[i]) & (lgms < lgms_bins[i+1])
        isred = is_red(lgms[sel], gcolor[sel])
        nred = float(np.sum(isred))
        print lgms_cens[i],
        # print np.sum(sel),
        print nred,
        isblue = ~isred
        nblue = float(np.sum(isblue))
        print nblue
        if nblue > 5:
            # plt.hist(gcolor[sel][isred], bins=30)
            # plt.hist(gcolor[sel][isblue], bins=30)
            # plt.show()
            pass
        lgms_cens_red[i] = np.log10(np.mean(10**lgms[sel][isred]))
        lgmh_cens_red[i] = np.mean(lgmh[sel][isred])
        lgmh_err_red[i] = np.std(lgmh[sel][isred])/np.sqrt(nred)
        lgms_cens_blue[i] = np.log10(np.mean(10**lgms[sel][isblue]))
        lgmh_cens_blue[i] = np.mean(lgmh[sel][isblue])
        lgmh_err_blue[i] = np.std(lgmh[sel][isblue])/np.sqrt(nblue)
    if has_fred:
        h = 0.701
        #
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
        # get HMF
        cosmo = CosmoParams(omega_M_0=0.27, sigma_8=0.82, h=0.70, omega_b_0=0.0469, n=0.95, set_flat=True)
        # Mh_arr, dndMh_arr = get_halofuncs(z=0.1, cosmo=cosmo, DELTA_HALO=200.0, mmin=1.e9, mmax=1.e16, nmbin=101)[:2]
        Mh_arr, dndMh_arr = read_mock_hmf(mockfile, mmin=1.e9, mmax=1.e16, nmbin=101, h=0.701)[:2]
        # get stellar mass grid
        # Ms_arr = np.logspace(8, 12, 201)
        Ms_arr = np.logspace(6, 15, 501)
        # normalize HMF
        # nhalo = np.trapz(dndMh_arr, x=Mh_arr)
        # p_Mh = dndMh_arr / nhalo
        # p_lnMh = p_Mh * Mh_arr
        #
        N_2darr = np.zeros((Ms_arr.size, Mh_arr.size))
        lnMs_mean_arr = shmr.log_stellarmass_mean(np.log(Mh_arr))
        sigma_arr = f_sigma_lnMs(Mh_arr)
        denom_arr = sigma_arr*np.sqrt(2.0*np.pi)
        lnMs_arr = np.log(Ms_arr)
        for i in xrange(Mh_arr.size):
            # at each log bin
            N_2darr[:, i] = np.exp(-0.5*((lnMs_arr - lnMs_mean_arr[i])/sigma_arr[i])**2) / denom_arr[i]
            # print np.trapz(N_2darr[:, i], x=lnMs_arr),
        #
        quenching = 'GD15'
        lgmhqc = 11.94779
        lgmhqs = 12.34138
        muc = 0.41160
        mus = 0.24031
        rfg = RedFractionGenerator(quenching=quenching, lgmhqc=lgmhqc, muc=muc,
                                   lgmhqs=lgmhqs, mus=mus)
        f_fred = rfg.make_fred_funcs(set_grid=True)[0]
        Nred_2darr = N_2darr * f_fred(Ms_arr, Mh_arr)
        Nblue_2darr = N_2darr * (1.0 - f_fred(Ms_arr, Mh_arr))
        hsmr_red = HSMR_split(Ms_arr, Mh_arr, Nred_2darr, dndMh_arr)
        hsmr_blue = HSMR_split(Ms_arr, Mh_arr, Nblue_2darr, dndMh_arr)
        # output
        lgmh_red = np.zeros(lgms_cens.size)
        lgms_red = np.zeros(lgms_cens.size)
        lgmh_blue = np.zeros(lgms_cens.size)
        lgms_blue = np.zeros(lgms_cens.size)
        for i in xrange(lgms_cens.size):
            Ms0 = 10**lgms_bins[i]
            Ms1 = 10**lgms_bins[i+1]
            # Ms_avg, _lnMh_mean, _lnMh_mean1, _lnMh_mean2
            lgms_red[i], lgmh_red[i] = hsmr_red.get_Mh_at_Msbin(Ms0, Ms1)[:2]
            lgms_red[i] = np.log10(lgms_red[i])
            lgmh_red[i] = lgmh_red[i] / np.log(10.0) + np.log10(h)
            lgms_blue[i], lgmh_blue[i] = hsmr_blue.get_Mh_at_Msbin(Ms0, Ms1)[:2]
            lgms_blue[i] = np.log10(lgms_blue[i])
            lgmh_blue[i] = lgmh_blue[i] / np.log(10.0) + np.log10(h)
        plt.plot(lgms_red, lgmh_red, 'r--')
        plt.plot(lgms_blue, lgmh_blue, 'b--')
        #
    # artificially increase the errorbar
    lgmh_err_red += 0.05
    lgmh_err_blue += 0.15
    plt.errorbar(lgms_cens_red, lgmh_cens_red,  yerr=lgmh_err_red, marker="o", ms=5, color="r")
    plt.errorbar(lgms_cens_blue, lgmh_cens_blue,  yerr=lgmh_err_blue, marker="s", ms=5, color="b")
    plt.xlabel(r"$\lg\;M_*\;[M_\odot/h^2]$")
    plt.ylabel(r"$\lg\;M_h\;[M_\odot/h]$")
    plt.xlim(9.8, 12)
    plt.show()

def test_mock_shmr_split(mockfile):
    """Check the color-split version of stellar to halo mass relation in the mock."""
    galrec = read_mock(mockfile)
    iscen = galrec['lg_halo_mass'] > 1
    lgmh = galrec['lg_halo_mass'][iscen]
    lgms = galrec['lg_stellar_mass'][iscen]
    gcolor = galrec['g-r'][iscen]
    lgmh_bins = np.linspace(11.4, 15.0, 35)
    lgmh_cens = (lgmh_bins[1:] + lgmh_bins[:-1]) / 2.0
    lgms_cens_red = np.empty_like(lgmh_cens)
    lgms_scas_red = np.empty_like(lgmh_cens)
    lgms_cens_blue = np.empty_like(lgmh_cens)
    lgms_scas_blue = np.empty_like(lgmh_cens)
    for i in xrange(lgmh_cens.size):
        sel = (lgmh >= lgmh_bins[i]) & (lgmh < lgmh_bins[i+1])
        isred = is_red(lgms[sel], gcolor[sel])
        lgms_cens_red[i] = np.mean(lgms[sel][isred])
        lgms_scas_red[i] = np.std(lgms[sel][isred])
        isblue = ~isred
        lgms_cens_blue[i] = np.mean(lgms[sel][isblue])
        lgms_scas_blue[i] = np.std(lgms[sel][isblue])
    if has_fred:
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
        plt.plot(lgmh_cens, lgmsarr, 'k-')
        plt.plot(lgmh_cens, lgmsarr + lgmssca , 'k--')
        plt.plot(lgmh_cens, lgmsarr - lgmssca , 'k--')
    plt.plot(lgmh_cens, lgms_cens_red, 'r-')
    plt.plot(lgmh_cens, lgms_cens_red+lgms_scas_red, 'r--')
    plt.plot(lgmh_cens, lgms_cens_red-lgms_scas_red, 'r--')
    plt.plot(lgmh_cens, lgms_cens_blue, 'b-')
    plt.plot(lgmh_cens, lgms_cens_blue+lgms_scas_blue, 'b--')
    plt.plot(lgmh_cens, lgms_cens_blue-lgms_scas_blue, 'b--')
    plt.xlabel(r"$M_h\;[M_\odot/h]$")
    plt.ylabel(r"$M_*\;[M_\odot/h^2]$")
    plt.show()

if __name__ == "__main__":
    # mockfile = '/Users/ying/Dropbox/Public/iHODcatalog_bolshoi.h5'
    mockfile = '/Users/ying/Data/ihodmock/standard/iHODcatalog_bolshoi.h5'
    # mockfile = '/Users/ying/Data/ihodmock/standard/iHODcatalog_mdr1.h5'
    # test_mock_shmr_split(mockfile)
    test_mock_hsmr_split(mockfile)

