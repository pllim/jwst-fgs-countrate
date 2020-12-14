import copy
import itertools

import numpy as np
import pandas as pd
import pytest

from fgscountrate.fgs_countrate_core import FGSCountrate
from fgscountrate import conversions


# Create data
# Aiming to match output from: fgs = FGSCountrate(guide_star_id="N13I000018"); data = fgs.query_gsc()
values = ['N13I000018', 420900912, 273.207, 65.5335, 8.30302e-05, 0.000185965,
          14.9447, 0.285722, 14.0877, 0.2927929, 13.7468, 0.239294,
          13.339, 0.0250000003, 12.993, 0.0270000007, 12.901, 0.0270000007,
          15.78594, 0.005142, 14.6547, 0.003211281, 14.27808, 0.003273380,
          14.1443, 0.003414216, 14.1067, 0.00433389]
index = ['hstID', 'gsc1ID', 'ra', 'dec', 'raErr', 'decErr',
         'JpgMag', 'JpgMagErr', 'FpgMag', 'FpgMagErr', 'NpgMag', 'NpgMagErr',
         'tmassJMag', 'tmassJMagErr', 'tmassHMag', 'tmassHMagErr', 'tmassKsMag', 'tmassKsMagErr',
         'SDSSuMag', 'SDSSuMagErr', 'SDSSgMag', 'SDSSgMagErr', 'SDSSrMag', 'SDSSrMagErr',
         'SDSSiMag', 'SDSSiMagErr', 'SDSSzMag', 'SDSSzMagErr']

BASE_DATA = pd.Series(values, index=index)


def test_convert_mag_to_jhk():
    """
    Using a different technique for choosing the conversion method and
    comparing that output to the technique used in calc_jhk_mag.
    This test is done by finding all possible combinations of the
    magnitudes that could be used for the conversion and checking
    that they have the same output with both techniques.

    """

    full_list = ['JpgMag', 'FpgMag', 'NpgMag', 'tmassJMag', 'tmassHMag',
                 'tmassKsMag', 'SDSSuMag', 'SDSSgMag', 'SDSSrMag',
                 'SDSSiMag', 'SDSSzMag']

    for L in range(0, len(full_list) + 1):
        for subset in itertools.combinations(full_list, L):
            if subset == ():
                continue

            # Recreate data
            fgs = FGSCountrate(guide_star_id="N13I000018", guider=1)
            fgs.gsc_series = BASE_DATA

            # Compute conversion
            delete_list = list(set(full_list) - set(subset))
            data2 = copy.copy(fgs.gsc_series)
            for i in delete_list:
                data2[i] = -999
            try:
                fgs.calc_jhk_mag(data=data2)
                error = False
            except ValueError:
                error = True

            # Compare output to here
            method_names = []
            for i in ['tmassJMag', 'tmassHMag', 'tmassKsMag']:
                if i in list(subset):
                    method_name_test = "convert_tmass_to_jhk"

                elif set(['SDSSgMag', 'SDSSzMag']).issubset(subset):
                    method_name_test = "convert_sdssgz_to_jhk"
                elif set(['SDSSgMag', 'SDSSiMag']).issubset(subset):
                    method_name_test = "convert_sdssgi_to_jhk"
                elif set(['SDSSiMag', 'SDSSzMag']).issubset(subset):
                    method_name_test = "convert_sdssiz_to_jhk"

                elif set(['JpgMag', 'NpgMag']).issubset(subset):
                    method_name_test = "convert_gsc2bjin_to_jhk"
                elif set(['FpgMag', 'NpgMag']).issubset(subset):
                    method_name_test = "convert_gsc2rfin_to_jhk"
                elif set(['JpgMag', 'FpgMag']).issubset(subset):
                    method_name_test = "convert_gsc2bjrf_to_jhk"

                else:
                    method_name_test = 'cannot_convert_to_jhk'

                method_names.append(method_name_test)

                if method_name_test != getattr(fgs, '{}_convert_method'.format(i[5].lower())):
                    if error is False:
                        print(subset)
                        print("    **", error, method_name_test, getattr(fgs, '{}_convert_method'.format(i[5].lower())))

                # If you could compute the JHK mags, check the same conversion method as used
                if error is False:
                    failure_message = 'For input {} and band {}: the test called {} while the calc_jhk_mag ' \
                                      'method called {}'.format(list(subset), i[5], method_name_test,
                                                                getattr(fgs, '{}_convert_method'.format(i[5].lower())))

                    assert method_name_test == getattr(fgs, '{}_convert_method'.format(i[5].lower())), failure_message

            # If you can't compute the JHK mags, check methods agree it's not possible
            if error is True:
                assert 'cannot_convert_to_jhk' in method_names


testdata = [
    (15, 15, 15, 25, 25, 25, 'convert_tmass_to_jhk', 'convert_tmass_to_jhk', 'convert_tmass_to_jhk'),
    (-999, -999, -999, 25, 15, 15, 'convert_sdssiz_to_jhk', 'convert_sdssiz_to_jhk', 'convert_sdssiz_to_jhk'),
]
ids = ['tmass conversion with faint sdss', 'sdss-zi conversion because of faint g-band']
@pytest.mark.parametrize("jmag, hmag, kmag, gmag, zmag, imag, convert_j, convert_h, convert_k", testdata, ids=ids)
def test_check_band_below_faint_limits_pass(jmag, hmag, kmag, gmag, zmag, imag, convert_j, convert_h, convert_k):
    """
    Test that the checking of faint bands will in certain cases have the code
    choose a conversion method that is not the "best" case conversion because
    of the existense of faint stars.
    """
    # Edit base data for specific test
    data = copy.copy(BASE_DATA)
    data['tmassJMag'] = jmag
    data['tmassHMag'] = hmag
    data['tmassKsMag'] = kmag
    data['SDSSgMag'] = gmag
    data['SDSSzMag'] = zmag
    data['SDSSiMag'] = imag
    data['JpgMag'] = -999
    data['FpgMag'] = -999
    data['NpgMag'] = -999

    fgs = FGSCountrate(guide_star_id="N13I000018", guider=1)
    fgs.calc_jhk_mag(data)

    assert fgs.j_convert_method == convert_j
    assert fgs.h_convert_method == convert_h
    assert fgs.k_convert_method == convert_k


def test_check_band_below_faint_limits_failure():
    """
    Check that if you have 1 catalog and 2+ bands below the faint limit,
    you cannot calculate the JHK bands and thus the magnitude. This should
    raise an error. This case is that you have only SDSS data, and 2 of the
    3 bands are below the faint limit.
    """
    # Edit base data for specific test
    data = copy.copy(BASE_DATA)
    data['tmassJMag'] = -999
    data['tmassHMag'] = -999
    data['tmassKsMag'] = -999
    data['SDSSgMag'] = 25
    data['SDSSzMag'] = 25
    data['SDSSiMag'] = 15
    data['JpgMag'] = -999
    data['FpgMag'] = -999
    data['NpgMag'] = -999

    fgs = FGSCountrate(guide_star_id="N13I000018", guider=1)

    with pytest.raises(Exception) as e_info:
        fgs.calc_jhk_mag(data)
    assert 'There is not enough information on this guide star' in str(e_info.value)
    assert 'tmassJMag' in str(e_info.value)


def test_tmass_to_jhk():
    """
    Check the _tmass_to_jhk method produces the expected result
    which is no change to the input values
    """
    # Create instance and get data
    fgs = FGSCountrate(guide_star_id="N13I000018", guider=2)
    fgs.gsc_series = BASE_DATA

    # Change data
    input_j = 10
    input_h = 11
    input_k = 12
    input_j_err = 0.10
    input_h_err = 0.11
    input_k_err = 0.12
    fgs.gsc_series['tmassJMag'] = input_j
    fgs.gsc_series['tmassHMag'] = input_h
    fgs.gsc_series['tmassKsMag'] = input_k
    fgs.gsc_series['tmassJMagErr'] = input_j_err
    fgs.gsc_series['tmassHMagErr'] = input_h_err
    fgs.gsc_series['tmassKsMagErr'] = input_k_err

    # Run method with series input
    j_ser, j_err_ser = conversions.convert_tmass_to_jhk(data=fgs.gsc_series, output_mag='J')
    h_ser, h_err_ser = conversions.convert_tmass_to_jhk(data=fgs.gsc_series, output_mag='H')
    k_ser, k_err_ser = conversions.convert_tmass_to_jhk(data=fgs.gsc_series, output_mag='K')

    # Run method with tuple input
    j_tup, j_err_tup = conversions.convert_tmass_to_jhk(data=(input_j, input_j_err), output_mag='J')
    h_tup, h_err_tup = conversions.convert_tmass_to_jhk(data=(input_h, input_h_err), output_mag='H')
    k_tup, k_err_tup = conversions.convert_tmass_to_jhk(data=(input_k, input_k_err), output_mag='K')

    # Check tuple and series input produces same result
    assert j_tup == j_ser
    assert h_tup == h_ser
    assert k_tup == k_ser

    # Check conversion function matches hand-check here
    assert j_ser == input_j
    assert h_ser == input_h
    assert k_ser == input_k

    # Check uncertainties
    assert np.isclose(j_err_ser, input_j_err, 1e-5)
    assert np.isclose(h_err_ser, input_h_err, 1e-5)
    assert np.isclose(k_err_ser, input_k_err, 1e-5)


def test_sdssgz_to_jhk():
    """Check the _sdssgz_to_jhk method produces the expected result """

    # Create instance and get data
    fgs = FGSCountrate(guide_star_id="N13I000018", guider=1)
    fgs.gsc_series = BASE_DATA

    # Change data
    input_g = 10
    input_z = 11
    input_g_err = 0.10
    input_z_err = 0.11
    fgs.gsc_series['SDSSgMag'] = input_g
    fgs.gsc_series['SDSSzMag'] = input_z
    fgs.gsc_series['SDSSgMagErr'] = input_g_err
    fgs.gsc_series['SDSSzMagErr'] = input_z_err

    # Run method with series input
    j_ser, j_err_ser = conversions.convert_sdssgz_to_jhk(data=fgs.gsc_series, output_mag='J')
    h_ser, h_err_ser = conversions.convert_sdssgz_to_jhk(data=fgs.gsc_series, output_mag='H')
    k_ser, k_err_ser = conversions.convert_sdssgz_to_jhk(data=fgs.gsc_series, output_mag='K')

    # Run method with tuple input
    data = (input_g, input_g_err, input_z, input_z_err)
    j_tup, j_err_tup = conversions.convert_sdssgz_to_jhk(data=data, output_mag='J')
    h_tup, h_err_tup = conversions.convert_sdssgz_to_jhk(data=data, output_mag='H')
    k_tup, k_err_tup = conversions.convert_sdssgz_to_jhk(data=data, output_mag='K')

    # Do calculation
    val = input_g - input_z
    output_j = input_g - 0.59 - 1.54*val + 0.20*val**2 - 0.04*val**3 + 0.002*val**4
    output_h = input_g - 0.77 - 1.78*val + 0.08*val**2 - 0.04*val**3 + 0.009*val**4
    output_k = input_g - 0.87 - 1.70*val + 0.01*val**2 - 0.07*val**3 + 0.001*val**4

    # Check tuple and series input produces same result
    assert j_tup == j_ser
    assert h_tup == h_ser
    assert k_tup == k_ser

    # Check conversion function matches hand-check here
    assert j_ser == output_j
    assert h_ser == output_h
    assert k_ser == output_k

    # Check uncertainties
    assert np.isclose(j_err_ser, 0.27265120293170503, 1e-5)
    assert np.isclose(h_err_ser, 0.28721836344780965, 1e-5)
    assert np.isclose(k_err_ser, 0.2946740530383989, 1e-5)


def test_sdssgi_to_jhk():
    """Check the _sdssgi_to_jhk method produces the expected result """

    # Create instance and get data
    fgs = FGSCountrate(guide_star_id="N13I000018", guider=2)
    fgs.gsc_series = BASE_DATA

    # Change data
    input_g = 10
    input_i = 11
    input_g_err = 0.10
    input_i_err = 0.11
    fgs.gsc_series['SDSSgMag'] = input_g
    fgs.gsc_series['SDSSiMag'] = input_i
    fgs.gsc_series['SDSSgMagErr'] = input_g_err
    fgs.gsc_series['SDSSiMagErr'] = input_i_err

    # Run method with series input
    j_ser, j_err_ser = conversions.convert_sdssgi_to_jhk(data=fgs.gsc_series, output_mag='J')
    h_ser, h_err_ser = conversions.convert_sdssgi_to_jhk(data=fgs.gsc_series, output_mag='H')
    k_ser, k_err_ser = conversions.convert_sdssgi_to_jhk(data=fgs.gsc_series, output_mag='K')

    # Run method with tuple input
    data = (input_g, input_g_err, input_i, input_i_err)
    j_tup, j_err_tup = conversions.convert_sdssgi_to_jhk(data=data, output_mag='J')
    h_tup, h_err_tup = conversions.convert_sdssgi_to_jhk(data=data, output_mag='H')
    k_tup, k_err_tup = conversions.convert_sdssgi_to_jhk(data=data, output_mag='K')

    # Do calculation
    val = input_g - input_i
    output_j = input_g - 0.411 - 2.260*val + 0.826*val**2 - 0.317*val**3 + 0.037*val**4
    output_h = input_g - 0.597 - 2.400*val + 0.450*val**2 - 0.078*val**3 + 0.00025*val**4
    output_k = input_g - 0.637 - 2.519*val + 0.568*val**2 - 0.151*val**3 + 0.013*val**4

    # Check tuple and series input produces same result
    assert j_tup == j_ser
    assert h_tup == h_ser
    assert k_tup == k_ser

    # Check conversion function matches hand-check here
    assert j_ser == output_j
    assert h_ser == output_h
    assert k_ser == output_k

    # Check uncertainties
    assert np.isclose(j_err_ser, 0.700281564042489, 1e-5)
    assert np.isclose(h_err_ser, 0.4895968466769927, 1e-5)


def test_sdssiz_to_jhk():
    """Check the _sdssiz_to_jhk method produces the expected result """

    # Create instance and get data
    fgs = FGSCountrate(guide_star_id="N13I000018", guider=1)
    fgs.gsc_series = BASE_DATA

    # Change data
    input_i = 10
    input_z = 11
    input_i_err = 0.10
    input_z_err = 0.11
    fgs.gsc_series['SDSSiMag'] = input_i
    fgs.gsc_series['SDSSzMag'] = input_z
    fgs.gsc_series['SDSSiMagErr'] = input_i_err
    fgs.gsc_series['SDSSzMagErr'] = input_z_err

    # Run method with series input
    j_ser, j_err_ser = conversions.convert_sdssiz_to_jhk(data=fgs.gsc_series, output_mag='J')
    h_ser, h_err_ser = conversions.convert_sdssiz_to_jhk(data=fgs.gsc_series, output_mag='H')
    k_ser, k_err_ser = conversions.convert_sdssiz_to_jhk(data=fgs.gsc_series, output_mag='K')

    # Run method with tuple input
    data = (input_i, input_i_err, input_z, input_z_err)
    j_tup, j_err_tup = conversions.convert_sdssiz_to_jhk(data=data, output_mag='J')
    h_tup, h_err_tup = conversions.convert_sdssiz_to_jhk(data=data, output_mag='H')
    k_tup, k_err_tup = conversions.convert_sdssiz_to_jhk(data=data, output_mag='K')

    # Do calculation
    val = input_i - input_z
    output_j = input_i - 0.794 - 2.839*val + 3.071*val**2 - 3.139*val**3 + 1.164*val**4
    output_h = input_i - 1.051 - 5.361*val + 8.398*val**2 - 7.240*val**3 + 2.111*val**4
    output_k = input_i - 1.127 - 5.379*val + 6.454*val**2 - 3.499*val**3 + 0.057*val**4

    # Check tuple and series input produces same result
    assert j_tup == j_ser
    assert h_tup == h_ser
    assert k_tup == k_ser

    # Check conversion function matches hand-check here
    assert j_ser == output_j
    assert h_ser == output_h
    assert k_ser == output_k

    # Check uncertainties
    assert np.isclose(j_err_ser, 3.4386988607490627, 1e-5)
    assert np.isclose(h_err_ser, 7.868447911817022, 1e-5)
    assert np.isclose(k_err_ser, 4.308133116754037, 1e-5)


def test_gsc2bjin_to_jhk():
    """Check the _gsc2bjin_to_jhk method produces the expected result """

    # Create instance and get data
    fgs = FGSCountrate(guide_star_id="N13I000018", guider=2)
    fgs.gsc_series = BASE_DATA

    # Change data
    input_bj = 10
    input_in = 11
    input_bj_err = 0.10
    input_in_err = 0.11
    fgs.gsc_series['JpgMag'] = input_bj
    fgs.gsc_series['NpgMag'] = input_in
    fgs.gsc_series['JpgMagErr'] = input_bj_err
    fgs.gsc_series['NpgMagErr'] = input_in_err

    # Run method with series input
    j_ser, j_err_ser = conversions.convert_gsc2bjin_to_jhk(data=fgs.gsc_series, output_mag='J')
    h_ser, h_err_ser = conversions.convert_gsc2bjin_to_jhk(data=fgs.gsc_series, output_mag='H')
    k_ser, k_err_ser = conversions.convert_gsc2bjin_to_jhk(data=fgs.gsc_series, output_mag='K')

    # Run method with tuple input
    data = (input_bj, input_bj_err, input_in, input_in_err)
    j_tup, j_err_tup = conversions.convert_gsc2bjin_to_jhk(data=data, output_mag='J')
    h_tup, h_err_tup = conversions.convert_gsc2bjin_to_jhk(data=data, output_mag='H')
    k_tup, k_err_tup = conversions.convert_gsc2bjin_to_jhk(data=data, output_mag='K')

    # Do calculation
    val = input_bj - input_in
    output_j = input_bj - 1.30*val - 0.15
    output_h = input_bj + 0.06*val**2 - 1.71*val - 0.10
    output_k = input_bj + 0.06*val**2 - 1.78*val - 0.11

    # Check tuple and series input produces same result
    assert j_tup == j_ser
    assert h_tup == h_ser
    assert k_tup == k_ser

    # Check conversion function matches hand-check here
    assert j_ser == output_j
    assert h_ser == output_h
    assert k_ser == output_k

    # Check uncertainties
    assert np.isclose(j_err_ser, 0.24527127838375157, 1e-5)
    assert np.isclose(h_err_ser, 0.3236128314452318, 1e-5)
    assert np.isclose(k_err_ser, 0.3649689782378769, 1e-5)


def test_gsc2rfin_to_jhk():
    """Check the _gsc2rfin_to_jhk method produces the expected result """

    # Create instance and get data
    fgs = FGSCountrate(guide_star_id="N13I000018", guider=1)
    fgs.gsc_series = BASE_DATA

    # Change data
    input_rf = 10
    input_in = 11
    input_rf_err = 0.10
    input_in_err = 0.11
    fgs.gsc_series['FpgMag'] = input_rf
    fgs.gsc_series['NpgMag'] = input_in
    fgs.gsc_series['FpgMagErr'] = input_rf_err
    fgs.gsc_series['NpgMagErr'] = input_in_err

    # Run method with series input
    j_ser, j_err_ser = conversions.convert_gsc2rfin_to_jhk(data=fgs.gsc_series, output_mag='J')
    h_ser, h_err_ser = conversions.convert_gsc2rfin_to_jhk(data=fgs.gsc_series, output_mag='H')
    k_ser, k_err_ser = conversions.convert_gsc2rfin_to_jhk(data=fgs.gsc_series, output_mag='K')

    # Run method with tuple input
    data = (input_rf, input_rf_err, input_in, input_in_err)
    j_tup, j_err_tup = conversions.convert_gsc2rfin_to_jhk(data=data, output_mag='J')
    h_tup, h_err_tup = conversions.convert_gsc2rfin_to_jhk(data=data, output_mag='H')
    k_tup, k_err_tup = conversions.convert_gsc2rfin_to_jhk(data=data, output_mag='K')

    # Do calculation
    val = input_rf - input_in
    output_j = input_rf + 0.01*val**2 - 1.56*val - 0.44
    output_h = input_rf + 0.25*val**2 - 2.17*val - 0.67
    output_k = input_rf + 0.28*val**2 - 2.35*val - 0.73

    # Check tuple and series input produces same result
    assert j_tup == j_ser
    assert h_tup == h_ser
    assert k_tup == k_ser

    # Check conversion function matches hand-check here
    assert j_ser == output_j
    assert h_ser == output_h
    assert k_ser == output_k

    # Check uncertainties
    assert np.isclose(j_err_ser, 0.30678481748776193, 1e-5)
    assert np.isclose(h_err_ser, 0.46706206827893615, 1e-5)
    assert np.isclose(k_err_ser, 0.5290933057070359, 1e-5)


def test_gsc2bjrf_to_jhk():
    """Check the _gsc2bjrf_to_jhk method produces the expected result """

    # Create instance and get data
    fgs = FGSCountrate(guide_star_id="N13I000018", guider=2)
    fgs.gsc_series = BASE_DATA

    # Change data
    input_bj = 10
    input_rf = 11
    input_bj_err = 0.10
    input_rf_err = 0.11
    fgs.gsc_series['JpgMag'] = input_bj
    fgs.gsc_series['FpgMag'] = input_rf
    fgs.gsc_series['JpgMagErr'] = input_bj_err
    fgs.gsc_series['FpgMagErr'] = input_rf_err

    # Run method with series input
    j_ser, j_err_ser = conversions.convert_gsc2bjrf_to_jhk(data=fgs.gsc_series, output_mag='J')
    h_ser, h_err_ser = conversions.convert_gsc2bjrf_to_jhk(data=fgs.gsc_series, output_mag='H')
    k_ser, k_err_ser = conversions.convert_gsc2bjrf_to_jhk(data=fgs.gsc_series, output_mag='K')

    # Run method with tuple input
    data = (input_bj, input_bj_err, input_rf, input_rf_err)
    j_tup, j_err_tup = conversions.convert_gsc2bjrf_to_jhk(data=data, output_mag='J')
    h_tup, h_err_tup = conversions.convert_gsc2bjrf_to_jhk(data=data, output_mag='H')
    k_tup, k_err_tup = conversions.convert_gsc2bjrf_to_jhk(data=data, output_mag='K')

    # Do calculation
    val = input_bj - input_rf
    output_j = input_bj - 0.39*val**2 - 0.96*val - 0.55
    output_h = input_bj - 0.24*val**2 - 1.66*val - 0.41
    output_k = input_bj - 0.26*val**2 - 1.70*val - 0.45

    # Check tuple and series input produces same result
    assert j_tup == j_ser
    assert h_tup == h_ser
    assert k_tup == k_ser

    # Check conversion function matches hand-check here
    assert j_ser == output_j
    assert h_ser == output_h
    assert k_ser == output_k

    # Check uncertainties
    assert np.isclose(j_err_ser, 0.3268272426848776, 1e-5)
    assert np.isclose(h_err_ser, 0.387910756252002, 1e-5)
    assert np.isclose(k_err_ser, 0.44586387576927555, 1e-5)
