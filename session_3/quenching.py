# 02 May 2017 16:39:54
import sys
sys.path.insert(1, '../session_2/')
sys.path.insert(1, '../session_1/')
try:
    from hod.fred import RedFractionGenerator
    from zypy.zycosmo import get_halofuncs, get_cosmo, CosmoParams
    from zypy.zyutil import getColors
    from hod.shmr import SHMR_Leauthaud, HSMR, HSMR_split, HSMR_grid
    from hod.predict import get_f_sigma_lnMs
    has_fred = True
except ImportError:
    has_fred = False
import numpy as np
import matplotlib.pyplot as plt
from read_mock import read_mock, read_mock_hmf
from read_m16 import read_m16_mass



"""Implement Halo Quenching."""


def is_red(lgms, gcolor):
    """Division of galaxies into red vs blue based on their relative position on
    the color-mass diagram."""
    cut = 0.8*(lgms/10.5)**0.6
    isred = gcolor >= cut
    return(isred)

def test_mock_hsmr_split(mockfile):
    """Check the color-split version of halo to stellar mass relations in the mock."""
    galrec = read_mock(mockfile)
    iscen = galrec['lg_halo_mass'] > 0
    # select mock central galaxies
    lgmh = galrec['lg_halo_mass'][iscen]
    lgms = galrec['lg_stellar_mass'][iscen]
    gcolor = galrec['g-r'][iscen]
    lgms_bins = np.linspace(8.8, 12.5, 20)
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
        lgms_cens_red[i] = np.log10(np.mean(10**lgms[sel][isred]))
        lgmh_cens_red[i] = np.mean(lgmh[sel][isred])
        lgmh_err_red[i] = np.std(lgmh[sel][isred])/np.sqrt(nred)
        #
        print nblue
        lgms_cens_blue[i] = np.log10(np.mean(10**lgms[sel][isblue]))
        lgmh_cens_blue[i] = np.mean(lgmh[sel][isblue])
        lgmh_err_blue[i] = np.std(lgmh[sel][isblue])/np.sqrt(nblue)
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
        # get HMF from theory
        # Mh_arr, dndMh_arr = get_halofuncs(z=0.1, cosmo=cosmo, DELTA_HALO=200.0, mmin=1.e9, mmax=1.e16, nmbin=201)[:2]
        # get HMF from simulation
        Mh_arr, dndlnMh_arr = read_mock_hmf(mockfile, mmin=1.e9, mmax=1.e16, nmbin=101, h=h)[:2]
        # get stellar mass vs. halo mass grid
        Ms_arr = np.logspace(6, 13, 301)
        #
        N_2darr = np.zeros((Ms_arr.size, Mh_arr.size))
        lnMs_mean_arr = shmr.log_stellarmass_mean(np.log(Mh_arr))
        sigma_arr = f_sigma_lnMs(Mh_arr)
        denom_arr = sigma_arr*np.sqrt(2.0*np.pi)
        lnMs_arr = np.log(Ms_arr)
        for i in xrange(Mh_arr.size):
            N_2darr[:, i] = np.exp(-0.5*((lnMs_arr - lnMs_mean_arr[i])/sigma_arr[i])**2) / denom_arr[i]
        if True:
            # get fred from theory
            quenching = 'GD15'
            lgmhqc = 11.94779
            lgmhqs = 12.34138
            muc = 0.41160
            mus = 0.24031
            # muc *= 1.2
            # lgmhqc *= 1.1
            rfg = RedFractionGenerator(quenching=quenching, lgmhqc=lgmhqc, muc=muc,
                                       lgmhqs=lgmhqs, mus=mus)
            f_fred = rfg.make_fred_funcs(set_grid=True)[0]
            fred_2darr =  f_fred(Ms_arr, Mh_arr)
            fblue_2darr = 1.0 - fred_2darr
        else:
            # get fred from mock
            fred = test_mock_fred(mockfile, mhmin=1.e9*h, mhmax=1.e16*h, nmhbin=101)[1]
            fred_2darr = fred.repeat(Ms_arr.size).reshape((Mh_arr.size, Ms_arr.size)).T
            fred_2darr[fred_2darr == 0.0] = 1.0
            fblue_2darr = 1.0 - fred_2darr
        Nred_2darr = N_2darr * fred_2darr
        Nblue_2darr = N_2darr * fblue_2darr
        hsmr_red = HSMR_grid(Ms_arr, Mh_arr, Nred_2darr, dndlnMh_arr)
        hsmr_blue = HSMR_grid(Ms_arr, Mh_arr, Nblue_2darr, dndlnMh_arr)
        # output
        lgmh_red = np.zeros(lgms_cens.size)
        lgms_red = np.zeros(lgms_cens.size)
        lgmh_blue = np.zeros(lgms_cens.size)
        lgms_blue = np.zeros(lgms_cens.size)
        for i in xrange(lgms_cens.size):
            Ms0 = 10**lgms_bins[i]
            Ms1 = 10**lgms_bins[i+1]
            lgms_red[i], lgmh_red[i] = hsmr_red.get_Mh_at_Msbin(Ms0, Ms1)[:2]
            lgms_red[i] = np.log10(lgms_red[i])
            lgmh_red[i] = lgmh_red[i] / np.log(10.0) + np.log10(h)
            lgms_blue[i], lgmh_blue[i] = hsmr_blue.get_Mh_at_Msbin(Ms0, Ms1)[:2]
            lgms_blue[i] = np.log10(lgms_blue[i])
            lgmh_blue[i] = lgmh_blue[i] / np.log(10.0) + np.log10(h)
        plt.plot(lgms_red, lgmh_red, 'r--', label="Red Pred")
        plt.plot(lgms_blue, lgmh_blue, 'b--', label="Blue Pred")
    plt.errorbar(lgms_cens_red, lgmh_cens_red,  yerr=lgmh_err_red, marker="o", ms=5, color="r", ls='None', label="Red Mock")
    plt.errorbar(lgms_cens_blue, lgmh_cens_blue,  yerr=lgmh_err_blue, marker="s", ms=5, color="b", ls='None', label="Blue Mock")
    if True:
        lgms, lgmh, lgmherr = read_m16_mass(True)
        plt.errorbar(lgms, lgmh, yerr=lgmherr, marker="o", ms=5, mfc="w", color="magenta", alpha=0.5, label="Red M16")
        lgms, lgmh, lgmherr = read_m16_mass(False)
        plt.errorbar(lgms, lgmh, yerr=lgmherr, marker="s", ms=5, mfc="w", color="cyan", alpha=0.5, label="Blue M16")
    plt.legend(loc=2)
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
    lgmh_bins = np.linspace(10.4, 15.0, 35)
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
    plt.xlabel(r"$\lg\;M_h\;[M_\odot/h]$")
    plt.ylabel(r"$f_{red}$")
    plt.show()
    return(lgmh_cens, f_red)

def test_misc(mockfile):
    """test intermedaite steps"""
    cosmo = CosmoParams(omega_M_0=0.27, sigma_8=0.82, h=0.70, omega_b_0=0.0469, n=0.95, set_flat=True)
    h = cosmo.h
    #
    galrec = read_mock(mockfile)
    iscen = galrec['lg_halo_mass'] > 0
    # convert to Msun
    lgmh = galrec['lg_halo_mass'][iscen] - np.log10(h)
    lgms = galrec['lg_stellar_mass'][iscen]
    gcolor = galrec['g-r'][iscen]
    #
    # get HMF
    # Mh_arr, _dndMh_arr = get_halofuncs(z=0.1, cosmo=cosmo, DELTA_HALO=200.0, mmin=1.e9, mmax=1.e15, nmbin=81)[:2]
    Mh_arr, dndlnMh_arr = read_mock_hmf(mockfile, mmin=1.e9, mmax=1.e15, nmbin=81, h=h)[:2]
    #
    dndlgMh_arr = dndlnMh_arr * np.log(10.0)
    lgMh_arr = np.log10(Mh_arr)
    dlgMh = np.log10(Mh_arr[1]/Mh_arr[0])
    # get stellar mass vs. halo mass grid
    Ms_arr = np.logspace(6, 13, 260)
    dlgMs = np.log10(Ms_arr[1]/Ms_arr[0])
    dlnMs = np.log(Ms_arr[1]/Ms_arr[0])
    lnMs_arr = np.log(Ms_arr)
    lgMs_arr = np.log10(Ms_arr)
    # make bins that centered on lnMs_arr
    lnMs_bin_arr = np.zeros(lnMs_arr.size + 1)
    lnMs_bin_arr[1:] = lnMs_arr + 0.5*dlnMs
    lnMs_bin_arr[0] = lnMs_arr[0] - 0.5*dlnMs
    lgMs_bin_arr = lnMs_bin_arr / np.log(10.0)
    #
    N_2darr = np.zeros((Ms_arr.size, Mh_arr.size))
    Nred_2darr = np.zeros((Ms_arr.size, Mh_arr.size))
    Nblue_2darr = np.zeros((Ms_arr.size, Mh_arr.size))
    pred_2darr = np.zeros((Ms_arr.size, Mh_arr.size))
    pblue_2darr = np.zeros((Ms_arr.size, Mh_arr.size))
    _N_2darr = np.zeros((Ms_arr.size, Mh_arr.size))
    _Nred_2darr = np.zeros((Ms_arr.size, Mh_arr.size))
    _Nblue_2darr = np.zeros((Ms_arr.size, Mh_arr.size))
    _pred_2darr = np.zeros((Ms_arr.size, Mh_arr.size))
    _pblue_2darr = np.zeros((Ms_arr.size, Mh_arr.size))
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
    #
    colors = getColors(Mh_arr.size)
    set_plot_grid = False
    for i in xrange(Mh_arr.size):
        if set_plot_grid:
            if np.mod(i, 4) == 0 and Mh_arr[i] > 2e11 and Mh_arr[i] < 1e13:
                print lgMh_arr[i]
            else:
                continue
        sel = (lgmh >= (lgMh_arr[i] - 0.5*dlgMh)) & (lgmh < (lgMh_arr[i] + 0.5*dlgMh))
        nhalo = np.sum(sel) * 1.0
        # print nhalo,
        # print dndlgMh_arr[i] * dlgMh * (250.0 / h)**3
        N_2darr[:, i] = np.exp(-0.5*((lnMs_arr - lnMs_mean_arr[i])/sigma_arr[i])**2) / denom_arr[i]
        #
        lgmhqc = 11.94779
        muc = 0.41160
        _mh = 10**(lgMh_arr[i] - lgmhqc)
        f_red = 1.0 - np.exp(-_mh**muc)
        # fred from theory
        Nred_2darr[:, i] = N_2darr[:, i] * nhalo * dlnMs * f_red
        pred_2darr[:, i] = N_2darr[:, i] * f_red
        Nblue_2darr[:, i] = N_2darr[:, i] * nhalo * dlnMs * (1.0 - f_red)
        pblue_2darr[:, i] = N_2darr[:, i] * (1.0 - f_red)
        #
        if set_plot_grid:
            # plt.plot(lgMs_arr, N_2darr[:, i] * nhalo * dlnMs, ls='-', color=colors[i], lw=2, alpha=0.5)
            # plt.plot(lgMs_arr, N_2darr[:, i] * nhalo * dlnMs * f_red, ls='-', color=colors[i], lw=3, alpha=0.5)
            plt.plot(lgMs_arr, N_2darr[:, i] * nhalo * dlnMs * (1.-f_red), ls='-', color=colors[i], lw=3, alpha=0.5)
        isred = is_red(lgms[sel], gcolor[sel])
        isblue = ~isred
        if np.sum(isred) >= 1:
            # fred from mock
            _N_2darr[:, i] = np.histogram(lgms[sel], bins=lgMs_bin_arr, normed=False)[0]
            _Nred_2darr[:, i] = np.histogram(lgms[sel][isred], bins=lgMs_bin_arr, normed=False)[0]
            _pred_2darr[:, i] = _Nred_2darr[:, i] / dlnMs / nhalo
            _Nblue_2darr[:, i] = np.histogram(lgms[sel][isblue], bins=lgMs_bin_arr, normed=False)[0]
            _pblue_2darr[:, i] = _Nblue_2darr[:, i] / dlnMs / nhalo
            if set_plot_grid:
                # plt.plot(lgMs_arr, _N_2darr[:, i], ls='-', color=colors[i], lw=1, alpha=0.8)
                # plt.plot(lgMs_arr, _Nred_2darr[:, i], ls='--', color=colors[i], lw=1, alpha=1.0)
                plt.plot(lgMs_arr, _Nblue_2darr[:, i], ls='--', color=colors[i], lw=1, alpha=1.0)
    if set_plot_grid:
        plt.xlabel(r'$M_*$')
        plt.yscale('log')
        plt.ylim(1, 4e3)
        plt.show()
        quit()
    #
    colors = getColors(Ms_arr.size)
    set_plot_dist = True
    #
    hsmr_blue = HSMR_grid(Ms_arr, Mh_arr, pblue_2darr, dndlnMh_arr)
    _hsmr_blue = HSMR_grid(Ms_arr, Mh_arr, _pblue_2darr, dndlnMh_arr)
    hsmr_red = HSMR_grid(Ms_arr, Mh_arr, pred_2darr, dndlnMh_arr)
    _hsmr_red = HSMR_grid(Ms_arr, Mh_arr, _pred_2darr, dndlnMh_arr)
    for j in xrange(Ms_arr.size):
        if set_plot_dist:
            if np.mod(j, 2) == 0 and Ms_arr[j] > 9e10 and Ms_arr[j] < 1e11:
                print np.log10(Ms_arr[j])
            else:
                continue
        else:
            # if Ms_arr[j] > 9.9e10 and Ms_arr[j] < 1e11:
            # if Ms_arr[j] > 5.8e10 and Ms_arr[j] < 6e10:
            # if Ms_arr[j] > 1.8e10 and Ms_arr[j] < 2e10:
            if Ms_arr[j] > 3.8e10 and Ms_arr[j] < 4e10:
            # if Ms_arr[j] > 4.8e10 and Ms_arr[j] < 5e10:
            # if Ms_arr[j] > 1e10 and Ms_arr[j] < 1.1e10:
                pass
            else:
                continue
        if set_plot_dist:
            plt.plot(lgMh_arr, Nblue_2darr[j, :], ls='-', color=colors[j], lw=3, alpha=0.5)
            plt.plot(lgMh_arr, _Nblue_2darr[j, :], ls='--', color=colors[j], lw=1, alpha=1.0)
        print '\n stellar mass',
        print lgMs_arr[j]
        print ('direct')
        print np.sum(Nblue_2darr[j, :] * lgMh_arr) / float(np.sum(Nblue_2darr[j, :])),
        # print np.sum(_Nblue_2darr[j, :] * lgMh_arr) / np.sum(_Nblue_2darr[j, :])
        # print np.sum(Nred_2darr[j, :] * lgMh_arr) / np.sum(Nred_2darr[j, :]),
        # print np.sum(_Nred_2darr[j, :] * lgMh_arr) / np.sum(_Nred_2darr[j, :])
        #
        print ('\n hsmr')
        # print hsmr_blue.get_Mh_at_Msbin(Ms_arr[j]*0.99, Ms_arr[j]*1.01)[1] / np.log(10.),
        # print _hsmr_blue.get_Mh_at_Msbin(Ms_arr[j]*0.99, Ms_arr[j]*1.01)[1] / np.log(10.)
        print hsmr_blue.get_Mh_at_Msbin(Ms_arr[j]/10**(0.5*dlgMs), Ms_arr[j]*10**(0.5*dlgMs))[1] / np.log(10.),
        # print _hsmr_blue.get_Mh_at_Msbin(Ms_arr[j]/10**(0.5*dlgMs), Ms_arr[j]*10**(0.5*dlgMs))[1] / np.log(10.)
        # print hsmr_red.get_Mh_at_Msbin(Ms_arr[j]/10**(0.5*dlgMs), Ms_arr[j]*10**(0.5*dlgMs))[1] / np.log(10.),
        # print _hsmr_red.get_Mh_at_Msbin(Ms_arr[j]/10**(0.5*dlgMs), Ms_arr[j]*10**(0.5*dlgMs))[1] / np.log(10.)
        # print np.log10(hsmr_blue.get_Mh_at_Msbin(Ms_arr[j]/10**(0.5*dlgMs), Ms_arr[j]*10**(0.5*dlgMs))[0]) ,
        # print np.log10(_hsmr_blue.get_Mh_at_Msbin(Ms_arr[j]/10**(0.5*dlgMs), Ms_arr[j]*10**(0.5*dlgMs))[0])
        print '\n'
        sel = (lgms >= (lgMs_arr[j] - 0.5*dlgMs)) & (lgms < (lgMs_arr[j] + 0.5*dlgMs))
        isred = is_red(lgms[sel], gcolor[sel])
        isblue = ~isred
        if False:
            # print np.sum(isblue)
            plt.hist(lgmh[sel][isblue], bins=lgMh_arr, normed=False)
            # plt.hist(lgmh[sel][isred], bins=lgMh_arr, normed=False)
            #
            p = hsmr_blue.get_plgMh_at_Msbin(Ms_arr[j]/10**(0.5*dlgMs), Ms_arr[j]*10**(0.5*dlgMs))*dlgMh*np.sum(isblue)
            # p = hsmr_red.get_plgMh_at_Msbin(Ms_arr[j]/10**(0.5*dlgMs), Ms_arr[j]*10**(0.5*dlgMs))*dlgMh*np.sum(isred)
            plt.plot(lgMh_arr, p, 'r-')
            plt.show()
    #
    #
    # print hsmr_blue.get_Mh_at_Msbin(8e10, 1e11)[1] / np.log(10.)
    # print _hsmr_blue.get_Mh_at_Msbin(8e10, 1e11)[1] / np.log(10.)
    if set_plot_dist:
        plt.xlabel(r'$M_h$')
        plt.yscale('log')
        plt.xlim(11, 15)
        plt.ylim(1, 2e4)
        plt.show()
        quit()


if __name__ == "__main__":
    # mockfile = '/Users/ying/Dropbox/Public/iHODcatalog_bolshoi.h5'
    mockfile = '/Users/ying/Data/ihodmock/standard/iHODcatalog_bolshoi.h5'
    # test_mock_shmr_split(mockfile)
    # test_mock_fred(mockfile)
    test_mock_hsmr_split(mockfile)
