# 02 May 2017 16:39:54
import sys
sys.path.insert(1, '../session_2/')
try:
    from hod.fred import RedFractionGenerator
    from zypy.zycosmo import get_halofuncs, get_cosmo, CosmoParams
    from zypy.zyutil import getColors
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

def is_blue(lgms, gcolor):
    cut = 0.8*(lgms/10.5)**0.6
    isblue = gcolor < cut
    return(isblue)

def test_mock_hsmr_split(mockfile):
    """Check the color-split version of halo to stellar mass relations in the mock."""
    galrec = read_mock(mockfile)
    iscen = galrec['lg_halo_mass'] > 0
    # select mock central galaxies
    lgmh = galrec['lg_halo_mass'][iscen]
    lgms = galrec['lg_stellar_mass'][iscen]
    gcolor = galrec['g-r'][iscen]
    lgms_bins = np.linspace(9.8, 12.5, 20)
    # lgms_bins = np.linspace(9.8, 12.5, 7)
    # lgms_bins = np.array([9.5, 9.8, 10.0, 10.2, 10.4, 10.6, 10.8, 11.0, 11.2, 11.4, 12.0])
    lgms_cens = (lgms_bins[1:] + lgms_bins[:-1]) / 2.0
    # initiate arrays
    # red
    lgms_cens_red = np.empty_like(lgms_cens)
    lgmh_cens_red = np.empty_like(lgms_cens)
    lgmh_err_red = np.empty_like(lgms_cens)
    # blue
    lgms_cens_blue = np.empty_like(lgms_cens)
    lgmh_cens_blue = np.empty_like(lgms_cens)
    lgmh_err_blue = np.empty_like(lgms_cens)
    # go thru each stellar mass bin
    for i in xrange(lgms_cens.size):
        print lgms_cens[i],
        sel = (lgms >= lgms_bins[i]) & (lgms < lgms_bins[i+1])
        isred = is_red(lgms[sel], gcolor[sel])
        nred = float(np.sum(isred))
        isblue = ~isred
        nblue = float(np.sum(isblue))
        print nred,
        print nblue
        lgms_cens_red[i] = np.log10(np.mean(10**lgms[sel][isred]))
        # lgms_cens_red[i] = np.mean(lgms[sel][isred])
        lgmh_cens_red[i] = np.mean(lgmh[sel][isred])
        # lgmh_cens_red[i] = np.median(lgmh[sel][isred])
        lgmh_err_red[i] = np.std(lgmh[sel][isred])/np.sqrt(nred)
        #
        lgms_cens_blue[i] = np.log10(np.mean(10**lgms[sel][isblue]))
        # lgms_cens_blue[i] = np.mean(lgms[sel][isblue])
        lgmh_cens_blue[i] = np.mean(lgmh[sel][isblue])
        # lgmh_cens_blue[i] = np.median(lgmh[sel][isblue])
        lgmh_err_blue[i] = np.std(lgmh[sel][isblue])/np.sqrt(nblue)
        if nblue > 2:
            # plt.hist(gcolor[sel][isred], bins=30)
            # plt.hist(gcolor[sel][isblue], bins=30)
            #
            plt.hist(lgmh[sel][isred], bins=30, color="red", alpha=0.5, normed=True)
            plt.hist(lgmh[sel][isblue], bins=30, color="blue", alpha=0.5, normed=True)
            plt.axvline(lgmh_cens_red[i], color="red")
            plt.axvline(lgmh_cens_blue[i], color="blue")
            plt.show()
            pass
    if has_fred:
        cosmo = CosmoParams(omega_M_0=0.27, sigma_8=0.82, h=0.70, omega_b_0=0.0469, n=0.95, set_flat=True)
        h = cosmo.h
        # shmr
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
        # Mh_arr, dndMh_arr = get_halofuncs(z=0.1, cosmo=cosmo, DELTA_HALO=200.0, mmin=1.e9, mmax=1.e16, nmbin=101)[:2]
        Mh_arr, dndMh_arr = read_mock_hmf(mockfile, mmin=1.e9, mmax=1.e16, nmbin=101, h=h)[:2]
        # get stellar mass vs. halo mass grid
        Ms_arr = np.logspace(9, 15, 301)
        #
        N_2darr = np.zeros((Ms_arr.size, Mh_arr.size))
        lnMs_mean_arr = shmr.log_stellarmass_mean(np.log(Mh_arr))
        sigma_arr = f_sigma_lnMs(Mh_arr)
        denom_arr = sigma_arr*np.sqrt(2.0*np.pi)
        lnMs_arr = np.log(Ms_arr)
        for i in xrange(Mh_arr.size):
            # at each log bin
            N_2darr[:, i] = np.exp(-0.5*((lnMs_arr - lnMs_mean_arr[i])/sigma_arr[i])**2) / denom_arr[i]
            jcut = np.searchsorted(Ms_arr, 1e10)
            # _norm =  np.trapz(N_2darr[:, i], x=lnMs_arr)
            N_2darr[:jcut, i] = 0
            # norm =  np.trapz(N_2darr[:, i], x=lnMs_arr)
            # N_2darr[:jcut, i] /= norm
        if False:
            quenching = 'GD15'
            lgmhqc = 11.94779
            lgmhqs = 12.34138
            muc = 0.41160
            mus = 0.24031
            rfg = RedFractionGenerator(quenching=quenching, lgmhqc=lgmhqc, muc=muc,
                                       lgmhqs=lgmhqs, mus=mus)
            f_fred = rfg.make_fred_funcs(set_grid=True)[0]
            fred_2darr =  f_fred(Ms_arr, Mh_arr)
            fblue_2darr = 1.0 - fred_2darr
        else:
            fred = test_mock_fred(mockfile, mhmin=1.e9*h, mhmax=1.e16*h, nmhbin=101)[1]
            fred_2darr = fred.repeat(Ms_arr.size).reshape((Mh_arr.size, Ms_arr.size)).T
            fred_2darr[fred_2darr == 0.0] = 1.0
            fblue_2darr = 1.0 - fred_2darr
        Nred_2darr = N_2darr * fred_2darr
        Nblue_2darr = N_2darr * fblue_2darr
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
        # plt.plot(np.log10(Ms_arr), hsmr_red.lnMh_med/np.log(10.0) + np.log10(h), 'r:')
        # plt.plot(np.log10(Ms_arr), hsmr_blue.lnMh_med/np.log(10.0) + np.log10(h), 'b:')
        # plt.plot(np.log10(Ms_arr), hsmr_red.lnMh_mean/np.log(10.0) + np.log10(h), 'r:')
        # plt.plot(np.log10(Ms_arr), hsmr_blue.lnMh_mean/np.log(10.0) + np.log10(h), 'b:')
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

def test_mock_fred(mockfile, mhmin=1e9, mhmax=1e15, nmhbin=35):
    """Check the red galaxy fraction in the mock."""
    galrec = read_mock(mockfile)
    iscen = galrec['lg_halo_mass'] > 1
    lgmh = galrec['lg_halo_mass'][iscen]
    lgms = galrec['lg_stellar_mass'][iscen]
    gcolor = galrec['g-r'][iscen]
    lgmhmin = np.log10(mhmin)
    lgmhmax = np.log10(mhmax)
    lgmh_bins = np.linspace(lgmhmin, lgmhmax, nmhbin)
    lgmh_cens = (lgmh_bins[1:] + lgmh_bins[:-1]) / 2.0
    f_red = np.zeros_like(lgmh_cens)
    for i in xrange(lgmh_cens.size):
        sel = (lgmh >= lgmh_bins[i]) & (lgmh < lgmh_bins[i+1])
        isred = is_red(lgms[sel], gcolor[sel])
        if np.sum(sel) > 1.0:
            f_red[i] = 1.0 * np.sum(isred) / (1. * np.sum(sel))
    if has_fred:
        h = 0.701
        lgmhqc = 11.94779
        muc = 0.41160
        _mh = 10**(lgmh_cens - np.log10(h) - lgmhqc)
        _f_red = 1.0 - np.exp(-_mh**muc)
        plt.plot(lgmh_cens, _f_red, 'k--', label="Theory")
        plt.axvline(lgmhqc - np.log10(h))
    plt.plot(lgmh_cens, f_red, 'r--', label="Mock")
    plt.legend(loc=2)
    # plt.plot(lgmh_cens, f_red/_f_red, 'r-', label="Mock")
    # plt.axhline(1)
    plt.xlabel(r"$\lg\;M_h\;[M_\odot/h]$")
    plt.ylabel(r"$f_{red}$")
    plt.show()
    return(lgmh_cens, f_red)

def test_grid(mockfile):
    """Check the color-split version of halo to stellar mass relations in the mock."""
    cosmo = CosmoParams(omega_M_0=0.27, sigma_8=0.82, h=0.70, omega_b_0=0.0469, n=0.95, set_flat=True)
    h = cosmo.h
    #
    galrec = read_mock(mockfile)
    iscen = galrec['lg_halo_mass'] > 0
    # select mock central galaxies
    lgmh = galrec['lg_halo_mass'][iscen] - np.log10(h)
    lgms = galrec['lg_stellar_mass'][iscen]
    gcolor = galrec['g-r'][iscen]
    #
    # get HMF
    # Mh_arr, dndMh_arr = get_halofuncs(z=0.1, cosmo=cosmo, DELTA_HALO=200.0, mmin=1.e9, mmax=1.e16, nmbin=101)[:2]
    Mh_arr, dndMh_arr = read_mock_hmf(mockfile, mmin=1.e11, mmax=1.e15, nmbin=8, h=h)[:2]
    lgMh_arr = np.log10(Mh_arr)
    dlgMh = np.log10(Mh_arr[1]/Mh_arr[0])
    # get stellar mass vs. halo mass grid
    Ms_arr = np.logspace(9, 13, 61)
    dlgMs = np.log10(Ms_arr[1]/Ms_arr[0])
    dlnMs = np.log(Ms_arr[1]/Ms_arr[0])
    #
    N_2darr = np.zeros((Ms_arr.size, Mh_arr.size))
    _N_2darr = np.zeros((Ms_arr.size, Mh_arr.size))
    #
    if True:
        # shmr
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
    lnMs_mean_arr = shmr.log_stellarmass_mean(np.log(Mh_arr))
    sigma_arr = f_sigma_lnMs(Mh_arr)
    denom_arr = sigma_arr*np.sqrt(2.0*np.pi)
    lnMs_arr = np.log(Ms_arr)
    lgMs_arr = np.log10(Ms_arr)
    lnMs_bin_arr = np.zeros(lnMs_arr.size + 1)
    lnMs_bin_arr[1:] = lnMs_arr + 0.5*dlnMs
    lnMs_bin_arr[0] = lnMs_arr[0] - 0.5*dlnMs
    lgMs_bin_arr = lnMs_bin_arr / np.log(10.0)
    colors = getColors(Mh_arr.size)
    for i in xrange(Mh_arr.size):
        sel = (lgmh >= (lgMh_arr[i] - 0.5*dlgMh)) & (lgmh < (lgMh_arr[i] + 0.5*dlgMh))
        nhalo = np.sum(sel) * 1.0
        # at each log bin
        print lgMh_arr[i],
        # if np.abs(lgMh_arr[i] - 11.2857142857) < 1e-1:
        if True:
            set_plot = True
        else:
            set_plot = False
        N_2darr[:, i] = np.exp(-0.5*((lnMs_arr - lnMs_mean_arr[i])/sigma_arr[i])**2) / denom_arr[i]
        print np.trapz(N_2darr[:, i], x=lnMs_arr),
        # plt.plot(lgMs_arr, N_2darr[:, i], ls='-', color=colors[i])
        if set_plot:
            # plt.plot(lgMs_arr, N_2darr[:, i] * dndMh_arr[i], ls='-', color=colors[i], lw=2, alpha=0.5)
            plt.plot(lgMs_arr, N_2darr[:, i] * nhalo, ls='-', color=colors[i], lw=2, alpha=0.5)
        # plt.plot(lgMs_arr, N_2darr[:, i] * dndMh_arr[i], 'r-', alpha=0.2)
        if np.sum(sel) >= 1:
            _N_2darr[:, i] = np.histogram(lgms[sel], bins=lgMs_bin_arr, normed=False)[0] / dlnMs / (np.sum(sel)*1.0)
            print np.trapz(_N_2darr[:, i], x=lnMs_arr)
            # plt.plot(lgMs_arr, _N_2darr[:, i], ls='--', color=colors[i])
            if set_plot:
                # plt.plot(lgMs_arr, _N_2darr[:, i] * dndMh_arr[i], ls='-', color=colors[i], lw=1, alpha=0.8)
                plt.plot(lgMs_arr, _N_2darr[:, i] * nhalo, ls='-', color=colors[i], lw=1, alpha=0.8)
            # plt.plot(lgMs_arr, _N_2darr[:, i] * dndMh_arr[i], 'k-', alpha=0.5)
        if set_plot:
            plt.xlabel(r'$M_*$')
            # plt.yscale('log')
            # plt.ylim(1e-30, 1e-10)
            # plt.ylim(1e-13, 1e1)
            plt.ylim(0, nhalo)
            plt.show()


if __name__ == "__main__":
    # mockfile = '/Users/ying/Dropbox/Public/iHODcatalog_bolshoi.h5'
    mockfile = '/Users/ying/Data/ihodmock/standard/iHODcatalog_bolshoi.h5'
    # mockfile = '/Users/ying/Data/ihodmock/standard/iHODcatalog_mdr1.h5'
    # test_mock_shmr_split(mockfile)
    # test_mock_fred(mockfile)
    # test_mock_hsmr_split(mockfile)
    test_grid(mockfile)

