# 01 May 2017 00:20:09
import os.path
import numpy as np
import matplotlib.pyplot as plt


"""Read Mandelbaum+2016 Weak Lensing data."""

m16path = '../data/M16/'



def read_m16(use_red=True, mass_bin='10.0_10.4'):
    """read data"""
    if use_red:
        # sm_10.0_10.4. - sm_10.4_10.7. - sm_10.7_11.0. - sm_11.0_11.2. - sm_11.2_11.4. - sm_11.4_11.6. - sm_11.6_15.0. - sm_11.0_15.0.
        fname = os.path.join(m16path, 'planck_lbg.ds.red.out')
        cols_dict ={
                '10.0_10.4': (0, 1, 2),
                '10.4_10.7': (0, 3, 4),
                '10.7_11.0': (0, 5, 6),
                '11.0_11.2': (0, 7, 8),
                '11.2_11.4': (0, 9, 10),
                '11.4_11.6': (0, 11, 12),
                '11.6_15.0': (0, 13, 14),
                '11.0_15.0': (0, 15, 16),
                }
    else:
        # sm_10.0_10.4. - sm_10.4_10.7. - sm_10.7_11.0. - sm_11.0_15.0.
        fname = os.path.join(m16path, 'planck_lbg.ds.blue.out')
        cols_dict ={
                '10.0_10.4': (0, 1, 2),
                '10.4_10.7': (0, 3, 4),
                '10.7_11.0': (0, 5, 6),
                '11.0_15.0': (0, 7, 8),
                }
    # Mpc/h, (h Msun/(physical pc)^2)
    rp, ds, ds_err = np.genfromtxt(fname, usecols=cols_dict[mass_bin], unpack=True)
    return(rp, ds, ds_err)


def test_read_m16(mass_bin="11.0_15.0"):
    rp, ds, ds_err = read_m16(use_red=True, mass_bin=mass_bin)
    plt.errorbar(rp, ds, yerr=ds_err, marker="o", ms=3, color="red")
    rp, ds, ds_err = read_m16(use_red=False, mass_bin=mass_bin)
    plt.errorbar(rp, ds, yerr=ds_err, marker="s", ms=3, color="blue")
    plt.xlabel(r"$R\;[Mpc/h]$")
    plt.ylabel(r"$\Delta\Sigma\;[h M_\odot/pc^2]$")
    plt.xscale('log')
    plt.yscale('log')
    plt.show()

if __name__ == "__main__":
    # test_read_m16(mass_bin='10.0_10.4')
    # test_read_m16(mass_bin='10.4_10.7')
    # test_read_m16(mass_bin='10.7_11.0')
    test_read_m16(mass_bin='11.0_15.0')

