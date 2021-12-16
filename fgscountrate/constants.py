"""Constants used in fgs_countrate_core.py and utils.py"""

# Universal constant
PLANCK = 6.625e-27

# Guide star catalog band information
GSC_BAND_NAMES = ['tmassJMag', 'tmassHMag', 'tmassKsMag',
                  'SDSSuMag', 'SDSSgMag', 'SDSSrMag', 'SDSSiMag',
                  'SDSSzMag', 'JpgMag', 'FpgMag', 'NpgMag']
GSC_BAND_WAVELENGTH = [1.25, 1.65, 2.17,
                       0.3551, 0.4680, 0.6166, 0.7480,
                       0.8932, 0.4660, 0.6450, 0.8500]

# Factor to use when calculating a band's missing uncertainty
BAND_ERR = 0.025

# FGS-Guider + OTE throughput
THROUGHPUT_G1 = {
    0.3551: 0.0,
    0.4660: 0.042,
    0.4680: 0.044,
    0.6166: 0.389,
    0.6450: 0.487,
    0.7480: 0.586,
    0.8500: 0.655,
    0.8932: 0.707,
    1.25: 0.688,
    1.65: 0.633,
    2.17: 0.723,
    3.0: 0.744,
    4.0: 0.690,
    5.0: 0.548,
    5.5: 0.041,
}
THROUGHPUT_G2 = {
    0.3551: 0.0,
    0.4660: 0.020,
    0.4680: 0.021,
    0.6166: 0.289,
    0.6450: 0.390,
    0.7480: 0.628,
    0.8500: 0.669,
    0.8932: 0.647,
    1.25: 0.761,
    1.65: 0.603,
    2.17: 0.635,
    3.0: 0.735,
    4.0: 0.738,
    5.0: 0.687,
    5.5: 0.040,
}

# Countrate conversion factors
CR_CONVERSION_G1 = 1.74
CR_CONVERSION_G2 = 1.57

# Magnitude conversion constant
MAG_CONVERSION_G1 = 28.29
MAG_CONVERSION_G2 = 28.20

# Conversion constants from Jmag to ABmag
ABMAG_CONSTANTS = {
    'tmassJMag': 0.90,
    'tmassHMag': 1.37,
    'tmassKsMag': 1.85,
    'SDSSuMag': 0.0,
    'SDSSgMag': 0.0,
    'SDSSrMag': 0.0,
    'SDSSiMag': 0.0,
    'SDSSzMag': 0.0,
    'JpgMag': -0.055,
    'FpgMag': 0.24,
    'NpgMag': 0.48,
}