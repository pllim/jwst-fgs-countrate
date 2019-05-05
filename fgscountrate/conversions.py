"""
Functions used to convert GSC, SDSS, and 2MASS magnitudes to J, H, and K mag.
"""

import pandas as pd


def convert_tmass_to_jhk(data, output_mag):
    """
    No conversion needed. 2MASS input is already in J,H,K band

    Parameters
    ----------
    data : Pandas Series or float
        Either a pandas Series as the output from a query
        of the Guide Star Catalog or a float of the
    output_mag : str
        The magnitude you want to convert to. Options
        are 'J', 'H', or 'K'.
    """

    if isinstance(data, pd.Series):
        j = data['tmassJmag']
        h = data['tmassHmag']
        k = data['tmassKsMag']
    elif isinstance(data, float):
        if output_mag.upper() == 'J':
            j = data
        elif output_mag.upper() == 'H':
            h = data
        elif output_mag.upper() == 'K':
            k = data
    else:
        raise TypeError("{} is not a valid type for data. Must be a float or a pd.Series output "
                        "from the Guide Star Catalog".format(type(data)))

    if output_mag.upper() == 'J':
        return j
    elif output_mag.upper() == 'H':
        return h
    elif output_mag.upper() == 'K':
        return k


def convert_sdssgz_to_jhk(data, output_mag):
    """
    Convert from SDSS_g mag and SDSS_z mag to J,H,K mag

    Parameters
    ----------
    data : Pandas Series or tuple
        Either a pandas Series as the output from a query
        of the Guide Star Catalog or a tuple containing the
        values (SDSS_g, SDSS_z)
    output_mag : str
        The magnitude you want to convert to. Options
        are 'J', 'H', or 'K'.
    """
    if isinstance(data, pd.Series):
        g = data['SDSSgMag']
        z = data['SDSSzMag']
    elif isinstance(data, tuple):
        g = data[0]
        z = data[1]
    else:
        raise TypeError("{} is not a valid type for data. Must be a tuple (g,z) or a pd.Series output "
                        "from the Guide Star Catalog".format(type(data)))

    if output_mag.upper() == 'J':
        j = g - 0.59 - 1.54 * (g - z) + 0.20 * (g - z) ** 2 - 0.04 * (g - z) ** 3 + 0.002 * (g - z) ** 4
        return j
    elif output_mag.upper() == 'H':
        h = g - 0.77 - 1.78 * (g - z) + 0.08 * (g - z) ** 2 - 0.04 * (g - z) ** 3 + 0.009 * (g - z) ** 4
        return h
    elif output_mag.upper() == 'K':
        k = g - 0.87 - 1.70 * (g - z) + 0.01 * (g - z) ** 2 - 0.07 * (g - z) ** 3 + 0.001 * (g - z) ** 4
        return k


def convert_sdssgi_to_jhk(data, output_mag):
    """
    Convert from SDSS_g mag and SDSS_i mag to J,H,K mag

    Parameters
    ----------
    data : Pandas Series or tuple
        Either a pandas Series as the output from a query
        of the Guide Star Catalog or a tuple containing the
        values (SDSS_g, SDSS_i)
    output_mag : str
        The magnitude you want to convert to. Options
        are 'J', 'H', or 'K'.
    """
    if isinstance(data, pd.Series):
        g = data['SDSSgMag']
        i = data['SDSSiMag']
    elif isinstance(data, tuple):
        g = data[0]
        i = data[1]
    else:
        raise TypeError("{} is not a valid type for data. Must be a tuple (g,i) or a pd.Series output "
                        "from the Guide Star Catalog".format(type(data)))

    if output_mag.upper() == 'J':
        j = g - 0.411 - 2.260 * (g - i) + 0.826 * (g - i) ** 2 - 0.317 * (g - i) ** 3 + 0.037 * (g - i) ** 4
        return j
    elif output_mag.upper() == 'H':
        h = g - 0.597 - 2.400 * (g - i) + 0.450 * (g - i) ** 2 - 0.078 * (g - i) ** 3 + 0.00025 * (g - i) ** 4
        return h
    elif output_mag.upper() == 'K':
        k = g - 0.637 - 2.519 * (g - i) + 0.568 * (g - i) ** 2 - 0.151 * (g - i) ** 3 + 0.013 * (g - i) ** 4
        return k


def convert_sdssiz_to_jhk(data, output_mag):
    """
    Convert from SDSS_i mag and SDSS_z mag to J,H,K mag

    Parameters
    ----------
    data : Pandas Series or tuple
        Either a pandas Series as the output from a query
        of the Guide Star Catalog or a tuple containing the
        values (SDSS_i, SDSS_z)
    output_mag : str
        The magnitude you want to convert to. Options
        are 'J', 'H', or 'K'.
    """
    if isinstance(data, pd.Series):
        i = data['SDSSiMag']
        z = data['SDSSzMag']
    elif isinstance(data, tuple):
        i = data[0]
        z = data[1]
    else:
        raise TypeError("{} is not a valid type for data. Must be a tuple (i,z) or a pd.Series output "
                        "from the Guide Star Catalog".format(type(data)))

    if output_mag.upper() == 'J':
        j = i - 0.794 - 2.839 * (i - z) + 3.071 * (i - z) ** 2 - 3.139 * (i - z) ** 3 + 1.164 * (i - z) ** 4
        return j
    elif output_mag.upper() == 'H':
        h = i - 1.051 - 5.361 * (i - z) + 8.398 * (i - z) ** 2 - 7.240 * (i - z) ** 3 + 2.111 * (i - z) ** 4
        return h
    elif output_mag.upper() == 'K':
        k = i - 1.127 - 5.379 * (i - z) + 6.454 * (i - z) ** 2 - 3.499 * (i - z) ** 3 + 0.057 * (i - z) ** 4
        return k


def convert_gsc2bjin_to_jhk(data, output_mag):
    """
    Convert from GSC2_B_J mag and GSC2_I_N mag to J,H,K mag

    Parameters
    ----------
    data : Pandas Series or tuple
        Either a pandas Series as the output from a query
        of the Guide Star Catalog or a tuple containing the
        values (GSC2_B_J, GSC2_I_N)
    output_mag : str
        The magnitude you want to convert to. Options
        are 'J', 'H', or 'K'.
    """
    if isinstance(data, pd.Series):
        b_j = data['JpgMag']
        i_n = data['NpgMag']
    elif isinstance(data, tuple):
        b_j = data[0]
        i_n = data[1]
    else:
        raise TypeError("{} is not a valid type for data. Must be a tuple (b_j,i_n) or a pd.Series output "
                        "from the Guide Star Catalog".format(type(data)))

    if output_mag.upper() == 'J':
        j = b_j - 1.30 * (b_j - i_n) - 0.15
        return j
    elif output_mag.upper() == 'H':
        h = b_j + 0.06 * (b_j - i_n) ** 2 - 1.71 * (b_j - i_n) - 0.10
        return h
    elif output_mag.upper() == 'K':
        k = b_j + 0.06 * (b_j - i_n) ** 2 - 1.78 * (b_j - i_n) - 0.11
        return k


def convert_gsc2rfin_to_jhk(data, output_mag):
    """
    Convert from GSC2_R_F mag and GSC2_I_N mag to J,H,K mag

    Parameters
    ----------
    data : Pandas Series or tuple
        Either a pandas Series as the output from a query
        of the Guide Star Catalog or a tuple containing the
        values (GSC2_R_F, GSC2_I_N)
    output_mag : str
        The magnitude you want to convert to. Options
        are 'J', 'H', or 'K'.
    """
    if isinstance(data, pd.Series):
        r_f = data['FpgMag']
        i_n = data['NpgMag']
    elif isinstance(data, tuple):
        r_f = data[0]
        i_n = data[1]
    else:
        raise TypeError("{} is not a valid type for data. Must be a tuple (r_f,i_n) or a pd.Series output "
                        "from the Guide Star Catalog".format(type(data)))

    if output_mag.upper() == 'J':
        j = r_f + 0.01 * (r_f - i_n) ** 2 - 1.56 * (r_f - i_n) - 0.44
        return j
    elif output_mag.upper() == 'H':
        h = r_f + 0.25 * (r_f - i_n) ** 2 - 2.17 * (r_f - i_n) - 0.67
        return h
    elif output_mag.upper() == 'K':
        k = r_f + 0.28 * (r_f - i_n) ** 2 - 2.35 * (r_f - i_n) - 0.73
        return k


def convert_gsc2bjrf_to_jhk(data, output_mag):
    """
    Convert from GSC2_B_J mag and GSC2_R_F mag to J,H,K mag

    Parameters
    ----------
    data : Pandas series or tuple
        Either a pandas Series as the output from a query
        of the Guide Star Catalog or a tuple containing the
        values (GSC2_B_J, GSC2_R_F)
    output_mag : str
        The magnitude you want to convert to. Options
        are 'J', 'H', or 'K'.
    """
    if isinstance(data, pd.Series):
        b_j = data['JpgMag']
        r_f = data['FpgMag']
    elif isinstance(data, tuple):
        b_j = data[0]
        r_f = data[1]
    else:
        raise TypeError("{} is not a valid type for data. Must be a tuple (b_j,r_f) or a pd.Series output "
                        "from the Guide Star Catalog".format(type(data)))

    if output_mag.upper() == 'J':
        j = b_j - 0.39 * (b_j - r_f) ** 2 - 0.96 * (b_j - r_f) - 0.55
        return j
    elif output_mag.upper() == 'H':
        h = b_j - 0.24 * (b_j - r_f) ** 2 - 1.66 * (b_j - r_f) - 0.41
        return h
    elif output_mag.upper() == 'K':
        k = b_j - 0.26 * (b_j - r_f) ** 2 - 1.70 * (b_j - r_f) - 0.45
        return k
