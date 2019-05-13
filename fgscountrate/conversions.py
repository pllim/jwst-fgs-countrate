"""
Functions used to convert GSC, SDSS, and 2MASS magnitudes to J, H, and K mag.

The sigma values representing the uncertainties in the transformations are
currently the robust standard deviations for the transformations.
"""

import numpy as np
import pandas as pd


def convert_tmass_to_jhk(data, output_mag):
    """
    Convert from 2MASS to J,H,K mag

    Note: No conversion is really needed here since the 2MASS
    input is already in J,H,K band. This function is really
    more of a place holder to be used in the the method
    FGSCountrate.convert_mag_to_jhk().

    Parameters
    ----------
    data : Pandas Series or tuple
        Either a pandas Series as the output from a query
        of the Guide Star Catalog or a a tuple containing the
        values (magnitude, error)
    output_mag : str
        The magnitude you want to convert to. Options
        are 'J', 'H', or 'K'.

    Returns
    -------
    JHK_mag : float
        The magnitude (J, H, or K based on the output_mag keyword)
    JHK_err : float
        The uncertainty of the magnitude
    """

    if isinstance(data, pd.Series):
        j_mag = data['tmassJmag']
        j_err = data['tmassJmagErr']
        h_mag = data['tmassHmag']
        h_err = data['tmassHmagErr']
        k_mag = data['tmassKsMag']
        k_err = data['tmassKsMagErr']
    elif isinstance(data, tuple):
        if output_mag.upper() == 'J':
            j_mag = data[0]
            j_err = data[1]
        elif output_mag.upper() == 'H':
            h_mag = data[0]
            h_err = data[1]
        elif output_mag.upper() == 'K':
            k_mag = data[0]
            k_err = data[1]
    else:
        raise TypeError("{} is not a valid type for data. Must be a tuple (mag, mag_err) or a pd.Series output "
                        "from the Guide Star Catalog".format(type(data)))

    if output_mag.upper() == 'J':
        return j_mag, j_err
    elif output_mag.upper() == 'H':
        return h_mag, h_err
    elif output_mag.upper() == 'K':
        return k_mag, k_err
    else:
        raise ValueError("output_mag must be set to either J, H, or K")


def convert_sdssgz_to_jhk(data, output_mag):
    """
    Convert from SDSS_g mag and SDSS_z mag to J,H,K mag

    Parameters
    ----------
    data : Pandas Series or tuple
        Either a pandas Series as the output from a query
        of the Guide Star Catalog or a tuple containing the
        values (SDSS_g, SDSS_g_err, SDSS_z, SDSS_z_err)
    output_mag : str
        The magnitude you want to convert to. Options
        are 'J', 'H', or 'K'.

    Returns
    -------
    JHK_mag : float
        The magnitude (J, H, or K based on the output_mag keyword)
    JHK_err : float
        The uncertainty of the magnitude
    """
    if isinstance(data, pd.Series):
        g_mag = data['SDSSgMag']
        g_err = data['SDSSgMagErr']
        z_mag = data['SDSSzMag']
        z_err = data['SDSSzMagErr']
    elif isinstance(data, tuple):
        g_mag = data[0]
        g_err = data[1]
        z_mag = data[2]
        z_err = data[3]
    else:
        raise TypeError("{} is not a valid type for data. Must be a tuple (g,g_err,z,z_err) or a pd.Series output "
                        "from the Guide Star Catalog".format(type(data)))

    if output_mag.upper() == 'J':
        def calc_j(g, z):
            return g - 0.59 - 1.54 * (g - z) + 0.20 * (g - z) ** 2 - 0.04 * (g - z) ** 3 + 0.002 * (g - z) ** 4
        j = calc_j(g_mag, z_mag)
        err_j_g = calc_j(g_mag + g_err, z_mag) - j
        err_j_z = calc_j(g_mag, z_mag + z_err) - j
        sigma_j_eqn = 0.1
        j_err = np.sqrt(err_j_g**2 + err_j_z**2 + sigma_j_eqn**2)
        return j, j_err
    elif output_mag.upper() == 'H':
        def calc_h(g, z):
            return g - 0.77 - 1.78 * (g - z) + 0.08 * (g - z) ** 2 - 0.04 * (g - z) ** 3 + 0.009 * (g - z) ** 4
        h = calc_h(g_mag, z_mag)
        err_h_g = calc_h(g_mag + g_err, z_mag) - h
        err_h_z = calc_h(g_mag, z_mag + z_err) - h
        sigma_h_eqn = 0.128
        h_err = np.sqrt(err_h_g**2 + err_h_z**2 + sigma_h_eqn**2)
        return h, h_err
    elif output_mag.upper() == 'K':
        def calc_k(g, z):
            return g - 0.87 - 1.70 * (g - z) + 0.01 * (g - z) ** 2 - 0.07 * (g - z) ** 3 + 0.001 * (g - z) ** 4
        k = calc_k(g_mag, z_mag)
        err_k_g = calc_k(g_mag + g_err, z_mag) - k
        err_k_z = calc_k(g_mag, z_mag + z_err) - k
        sigma_k_eqn = 0.179
        k_err = np.sqrt(err_k_g**2 + err_k_z**2 + sigma_k_eqn**2)
        return k, k_err
    else:
        raise ValueError("output_mag must be set to either J, H, or K")


def convert_sdssgi_to_jhk(data, output_mag):
    """
    Convert from SDSS_g mag and SDSS_i mag to J,H,K mag

    Parameters
    ----------
    data : Pandas Series or tuple
        Either a pandas Series as the output from a query
        of the Guide Star Catalog or a tuple containing the
        values (SDSS_g, SDSS_g_err, SDSS_i, SDSS_i_err)
    output_mag : str
        The magnitude you want to convert to. Options
        are 'J', 'H', or 'K'.

    Returns
    -------
    JHK_mag : float
        The magnitude (J, H, or K based on the output_mag keyword)
    JHK_err : float
        The uncertainty of the magnitude
    """
    if isinstance(data, pd.Series):
        g_mag = data['SDSSgMag']
        g_err = data['SDSSgMagErr']
        i_mag = data['SDSSiMag']
        i_err = data['SDSSiMagErr']
    elif isinstance(data, tuple):
        g_mag = data[0]
        g_err = data[1]
        i_mag = data[2]
        i_err = data[3]
    else:
        raise TypeError("{} is not a valid type for data. Must be a tuple (g,g_err,i,i_err) or a pd.Series output "
                        "from the Guide Star Catalog".format(type(data)))

    if output_mag.upper() == 'J':
        def calc_j(g, i):
            return g - 0.411 - 2.260 * (g - i) + 0.826 * (g - i) ** 2 - 0.317 * (g - i) ** 3 + 0.037 * (g - i) ** 4
        j = calc_j(g_mag, i_mag)
        err_j_g = calc_j(g_mag + g_err, i_mag) - j
        err_j_i = calc_j(g_mag, i_mag + i_err) - j
        sigma_j_eqn = 0.114
        j_err = np.sqrt(err_j_g**2 + err_j_i**2 + sigma_j_eqn**2)
        return j, j_err
    elif output_mag.upper() == 'H':
        def calc_h(g, i):
            return g - 0.597 - 2.400 * (g - i) + 0.450 * (g - i) ** 2 - 0.078 * (g - i) ** 3 + 0.00025 * (g - i) ** 4
        h = calc_h(g_mag, i_mag)
        err_h_g = calc_h(g_mag + g_err, i_mag) - h
        err_h_i = calc_h(g_mag, i_mag + i_err) - h
        sigma_h_eqn = 0.145
        h_err = np.sqrt(err_h_g**2 + err_h_i**2 + sigma_h_eqn**2)
        return h, h_err
    elif output_mag.upper() == 'K':
        def calc_k(g, i):
            return g - 0.637 - 2.519 * (g - i) + 0.568 * (g - i) ** 2 - 0.151 * (g - i) ** 3 + 0.013 * (g - i) ** 4
        k = calc_k(g_mag, i_mag)
        err_k_g = calc_k(g_mag + g_err, i_mag) - k
        err_k_i = calc_k(g_mag, i_mag + i_err) - k
        sigma_k_eqn = 0.203
        k_err = np.sqrt(err_k_g**2 + err_k_i**2 + sigma_k_eqn**2)
        return k, k_err
    else:
        raise ValueError("output_mag must be set to either J, H, or K")


def convert_sdssiz_to_jhk(data, output_mag):
    """
    Convert from SDSS_i mag and SDSS_z mag to J,H,K mag

    Parameters
    ----------
    data : Pandas Series or tuple
        Either a pandas Series as the output from a query
        of the Guide Star Catalog or a tuple containing the
        values (SDSS_i, SDSS_i_err, SDSS_z, SDSS_z_err)
    output_mag : str
        The magnitude you want to convert to. Options
        are 'J', 'H', or 'K'.

    Returns
    -------
    JHK_mag : float
        The magnitude (J, H, or K based on the output_mag keyword)
    JHK_err : float
        The uncertainty of the magnitude
    """
    if isinstance(data, pd.Series):
        i_mag = data['SDSSiMag']
        i_err = data['SDSSiMagErr']
        z_mag = data['SDSSzMag']
        z_err = data['SDSSzMagErr']
    elif isinstance(data, tuple):
        i_mag = data[0]
        i_err = data[1]
        z_mag = data[2]
        z_err = data[3]
    else:
        raise TypeError("{} is not a valid type for data. Must be a tuple (i,i_err,z,z_err) or a pd.Series output "
                        "from the Guide Star Catalog".format(type(data)))

    if output_mag.upper() == 'J':
        def calc_j(i, z):
            return i - 0.794 - 2.839 * (i - z) + 3.071 * (i - z) ** 2 - 3.139 * (i - z) ** 3 + 1.164 * (i - z) ** 4
        j = calc_j(i_mag, z_mag)
        err_j_i = calc_j(i_mag + i_err, z_mag) - j
        err_j_z = calc_j(i_mag, z_mag + z_err) - j
        sigma_j_eqn = 0.104
        j_err = np.sqrt(err_j_i**2 + err_j_z**2 + sigma_j_eqn**2)
        return j, j_err
    elif output_mag.upper() == 'H':
        def calc_h(i, z):
            return i - 1.051 - 5.361 * (i - z) + 8.398 * (i - z) ** 2 - 7.240 * (i - z) ** 3 + 2.111 * (i - z) ** 4
        h = calc_h(i_mag, z_mag)
        err_h_i = calc_h(i_mag + i_err, z_mag) - h
        err_h_z = calc_h(i_mag, z_mag + z_err) - h
        sigma_h_eqn = 0.141
        h_err = np.sqrt(err_h_i**2 + err_h_z**2 + sigma_h_eqn**2)
        return h, h_err
    elif output_mag.upper() == 'K':
        def calc_k(i, z):
            return i - 1.127 - 5.379 * (i - z) + 6.454 * (i - z) ** 2 - 3.499 * (i - z) ** 3 + 0.057 * (i - z) ** 4
        k = calc_k(i_mag, z_mag)
        err_k_i = calc_k(i_mag + i_err, z_mag) - k
        err_k_z = calc_k(i_mag, z_mag + z_err) - k
        sigma_k_eqn = 0.187
        k_err = np.sqrt(err_k_i**2 + err_k_z**2 + sigma_k_eqn**2)
        return k, k_err
    else:
        raise ValueError("output_mag must be set to either J, H, or K")


def convert_gsc2bjin_to_jhk(data, output_mag):
    """
    Convert from GSC2_B_J mag and GSC2_I_N mag to J,H,K mag

    Parameters
    ----------
    data : Pandas Series or tuple
        Either a pandas Series as the output from a query
        of the Guide Star Catalog or a tuple containing the
        values (GSC2_B_J, GSC2_B_J_err, GSC2_I_N, GSC2_I_N_err)
    output_mag : str
        The magnitude you want to convert to. Options
        are 'J', 'H', or 'K'.

    Returns
    -------
    JHK_mag : float
        The magnitude (J, H, or K based on the output_mag keyword)
    JHK_err : float
        The uncertainty of the magnitude
    """
    if isinstance(data, pd.Series):
        b_j_mag = data['JpgMag']
        b_j_err = data['JpgMagErr']
        i_n_mag = data['NpgMag']
        i_n_err = data['NpgMagErr']
    elif isinstance(data, tuple):
        b_j_mag = data[0]
        b_j_err = data[1]
        i_n_mag = data[2]
        i_n_err = data[3]
    else:
        raise TypeError("{} is not a valid type for data. Must be a tuple (b_j,b_j_err,i_n,i_n_err) "
                        "or a pd.Series output from the Guide Star Catalog".format(type(data)))

    if output_mag.upper() == 'J':
        def calc_j(b_j, i_n): return b_j - 1.30 * (b_j - i_n) - 0.15
        j = calc_j(b_j_mag, i_n_mag)
        err_j_b_j = calc_j(b_j_mag + b_j_err, i_n_mag) - j
        err_j_i_n = calc_j(b_j_mag, i_n_mag + i_n_err) - j
        sigma_j_eqn = 0.197
        j_err = np.sqrt(err_j_b_j**2 + err_j_i_n**2 + sigma_j_eqn**2)
        return j, j_err
    elif output_mag.upper() == 'H':
        def calc_h(b_j, i_n): return b_j + 0.06 * (b_j - i_n) ** 2 - 1.71 * (b_j - i_n) - 0.10
        h = calc_h(b_j_mag, i_n_mag)
        err_h_b_j = calc_h(b_j_mag + b_j_err, i_n_mag) - h
        err_h_i_n = calc_h(b_j_mag, i_n_mag + i_n_err) - h
        sigma_h_eqn = 0.239
        h_err = np.sqrt(err_h_b_j**2 + err_h_i_n**2 + sigma_h_eqn**2)
        return h, h_err
    elif output_mag.upper() == 'K':
        def calc_k(b_j, i_n): return b_j + 0.06 * (b_j - i_n) ** 2 - 1.78 * (b_j - i_n) - 0.11
        k = calc_k(b_j_mag, i_n_mag)
        err_k_b_j = calc_k(b_j_mag + b_j_err, i_n_mag) - k
        err_k_i_n = calc_k(b_j_mag, i_n_mag + i_n_err) - k
        sigma_k_eqn = 0.285
        k_err = np.sqrt(err_k_b_j**2 + err_k_i_n**2 + sigma_k_eqn**2)
        return k, k_err
    else:
        raise ValueError("output_mag must be set to either J, H, or K")


def convert_gsc2rfin_to_jhk(data, output_mag):
    """
    Convert from GSC2_R_F mag and GSC2_I_N mag to J,H,K mag

    Parameters
    ----------
    data : Pandas Series or tuple
        Either a pandas Series as the output from a query
        of the Guide Star Catalog or a tuple containing the
        values (GSC2_R_F, GSC2_R_F_err, GSC2_I_N, GSC2_I_N_err)
    output_mag : str
        The magnitude you want to convert to. Options
        are 'J', 'H', or 'K'.

    Returns
    -------
    JHK_mag : float
        The magnitude (J, H, or K based on the output_mag keyword)
    JHK_err : float
        The uncertainty of the magnitude
    """
    if isinstance(data, pd.Series):
        r_f_mag = data['FpgMag']
        r_f_err = data['FpgMagErr']
        i_n_mag = data['NpgMag']
        i_n_err = data['NpgMagErr']
    elif isinstance(data, tuple):
        r_f_mag = data[0]
        r_f_err = data[1]
        i_n_mag = data[2]
        i_n_err = data[3]
    else:
        raise TypeError("{} is not a valid type for data. Must be a tuple (r_f,r_f_err,i_n,i_n_err) or a "
                        "pd.Series output from the Guide Star Catalog".format(type(data)))

    if output_mag.upper() == 'J':
        def calc_j(r_f, i_n): return r_f + 0.01 * (r_f - i_n) ** 2 - 1.56 * (r_f - i_n) - 0.44
        j = calc_j(r_f_mag, i_n_mag)
        err_j_r_f = calc_j(r_f_mag + r_f_err, i_n_mag) - j
        err_j_i_n = calc_j(r_f_mag, i_n_mag + i_n_err) - j
        sigma_j_eqn = 0.246
        j_err = np.sqrt(err_j_r_f**2 + err_j_i_n**2 + sigma_j_eqn**2)
        return j, j_err
    elif output_mag.upper() == 'H':
        def calc_h(r_f, i_n): return r_f + 0.25 * (r_f - i_n) ** 2 - 2.17 * (r_f - i_n) - 0.67
        h = calc_h(r_f_mag, i_n_mag)
        err_h_r_f = calc_h(r_f_mag + r_f_err, i_n_mag) - h
        err_h_i_n = calc_h(r_f_mag, i_n_mag + i_n_err) - h
        sigma_h_eqn = 0.321
        h_err = np.sqrt(err_h_r_f**2 + err_h_i_n**2 + sigma_h_eqn**2)
        return h, h_err
    elif output_mag.upper() == 'K':
        def calc_k(r_f, i_n): return r_f + 0.28 * (r_f - i_n) ** 2 - 2.35 * (r_f - i_n) - 0.73
        k = calc_k(r_f_mag, i_n_mag)
        err_k_r_f = calc_k(r_f_mag + r_f_err, i_n_mag) - k
        err_k_i_n = calc_k(r_f_mag, i_n_mag + i_n_err) - k
        sigma_k_eqn = 0.374
        k_err = np.sqrt(err_k_r_f**2 + err_k_i_n**2 + sigma_k_eqn**2)
        return k, k_err
    else:
        raise ValueError("output_mag must be set to either J, H, or K")


def convert_gsc2bjrf_to_jhk(data, output_mag):
    """
    Convert from GSC2_B_J mag and GSC2_R_F mag to J,H,K mag

    Parameters
    ----------
    data : Pandas series or tuple
        Either a pandas Series as the output from a query
        of the Guide Star Catalog or a tuple containing the
        values (GSC2_B_J, GSC2_B_J_err, GSC2_R_F, GSC2_R_F_err)
    output_mag : str
        The magnitude you want to convert to. Options
        are 'J', 'H', or 'K'.

    Returns
    -------
    JHK_mag : float
        The magnitude (J, H, or K based on the output_mag keyword)
    JHK_err : float
        The uncertainty of the magnitude
    """
    if isinstance(data, pd.Series):
        b_j_mag = data['JpgMag']
        b_j_err = data['JpgMagErr']
        r_f_mag = data['FpgMag']
        r_f_err = data['FpgMagErr']
    elif isinstance(data, tuple):
        b_j_mag = data[0]
        b_j_err = data[1]
        r_f_mag = data[2]
        r_f_err = data[3]
    else:
        raise TypeError("{} is not a valid type for data. Must be a tuple (b_j,b_j_err,r_f,r_f_err) or a "
                        "pd.Series output from the Guide Star Catalog".format(type(data)))

    if output_mag.upper() == 'J':
        def calc_j(b_j, r_f): return b_j - 0.39 * (b_j - r_f) ** 2 - 0.96 * (b_j - r_f) - 0.55
        j = calc_j(b_j_mag, r_f_mag)
        err_j_b_j = calc_j(b_j_mag + b_j_err, r_f_mag) - j
        err_j_r_f = calc_j(b_j_mag, r_f_mag + r_f_err) - j
        sigma_j_eqn = 0.317
        j_err = np.sqrt(err_j_b_j**2 + err_j_r_f**2 + sigma_j_eqn**2)
        return j, j_err
    elif output_mag.upper() == 'H':
        def calc_h(b_j, r_f): return b_j - 0.24 * (b_j - r_f) ** 2 - 1.66 * (b_j - r_f) - 0.41
        h = calc_h(b_j_mag, r_f_mag)
        err_h_b_j = calc_h(b_j_mag + b_j_err, r_f_mag) - h
        err_h_r_f = calc_h(b_j_mag, r_f_mag + r_f_err) - h
        sigma_h_eqn = 0.366
        h_err = np.sqrt(err_h_b_j**2 + err_h_r_f**2 + sigma_h_eqn**2)
        return h, h_err
    elif output_mag.upper() == 'K':
        def calc_k(b_j, r_f): return b_j - 0.26 * (b_j - r_f) ** 2 - 1.70 * (b_j - r_f) - 0.45
        k = calc_k(b_j_mag, r_f_mag)
        err_k_b_j = calc_k(b_j_mag + b_j_err, r_f_mag) - k
        err_k_r_f = calc_k(b_j_mag, r_f_mag + r_f_err) - k
        sigma_k_eqn = 0.427
        k_err = np.sqrt(err_k_b_j**2 + err_k_r_f**2 + sigma_k_eqn**2)
        return k, k_err
    else:
        raise ValueError("output_mag must be set to either J, H, or K")
