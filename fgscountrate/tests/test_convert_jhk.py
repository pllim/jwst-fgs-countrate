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
            fgs.calc_jhk_mag(data=data2)

            # Compare output to here
            method_names = []
            for i in ['tmassJMag', 'tmassHMag', 'tmassKsMag']:
                if i in list(subset):
                    method_name_test = "convert_tmass_to_jhk"

                elif {'SDSSgMag', 'SDSSzMag'}.issubset(subset):
                    method_name_test = "convert_sdssgz_to_jhk"
                elif {'SDSSrMag', 'SDSSzMag'}.issubset(subset):
                    method_name_test = "convert_sdssrz_to_jhk"
                elif {'SDSSgMag', 'SDSSiMag'}.issubset(subset):
                    method_name_test = "convert_sdssgi_to_jhk"
                elif {'SDSSiMag', 'SDSSzMag'}.issubset(subset):
                    method_name_test = "convert_sdssiz_to_jhk"

                elif {'JpgMag', 'NpgMag'}.issubset(subset):
                    method_name_test = "convert_gsc2bjin_to_jhk"
                elif {'FpgMag', 'NpgMag'}.issubset(subset):
                    method_name_test = "convert_gsc2rfin_to_jhk"
                elif {'JpgMag', 'FpgMag'}.issubset(subset):
                    method_name_test = "convert_gsc2bjrf_to_jhk"

                else:
                    method_name_test = 'cannot_calculate'

                method_names.append(method_name_test)

                failure_message = f'For input {list(subset)} and band {i[5]}: the test called {method_name_test} ' \
                                  f'while the calc_jhk_mag method called ' \
                                  f'{getattr(fgs, f"{i[5].lower()}_convert_method")}'

                assert method_name_test == getattr(fgs, f'{i[5].lower()}_convert_method'), failure_message


testdata = [
    (15, 15, 15, 25, 25, 25, 25, 'convert_tmass_to_jhk', 'convert_tmass_to_jhk', 'convert_tmass_to_jhk'),
    (-999, -999, -999, 25, 15, 25, 15, 'convert_sdssiz_to_jhk', 'convert_sdssiz_to_jhk', 'convert_sdssiz_to_jhk'),
]
ids = ['tmass conversion with faint sdss', 'sdss-zi conversion because of faint g-band']
@pytest.mark.parametrize("jmag, hmag, kmag, gmag, zmag, rmag, imag, convert_j, convert_h, convert_k", testdata, ids=ids)
def test_check_band_below_faint_limits_pass(jmag, hmag, kmag, gmag, zmag, rmag, imag, convert_j, convert_h, convert_k):
    """
    Test that the checking of faint bands will in certain cases have the code
    choose a conversion method that is not the "best" case conversion because
    of the existence of faint stars.
    """
    # Edit base data for specific test
    data = copy.copy(BASE_DATA)
    data['tmassJMag'] = jmag
    data['tmassHMag'] = hmag
    data['tmassKsMag'] = kmag
    data['SDSSgMag'] = gmag
    data['SDSSrMag'] = rmag
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
    data['SDSSrMag'] = -999
    data['SDSSzMag'] = 25
    data['SDSSiMag'] = 15
    data['JpgMag'] = -999
    data['FpgMag'] = -999
    data['NpgMag'] = -999

    fgs = FGSCountrate(guide_star_id="N13I000018", guider=1)

    _ = fgs.calc_jhk_mag(data)

    assert fgs.j_convert_method == 'cannot_calculate'
    assert fgs.h_convert_method == 'cannot_calculate'
    assert fgs.k_convert_method == 'cannot_calculate'

    with pytest.raises(Exception) as e_info:
        fgs.calc_fgs_cr_mag_and_err()
    assert 'Cannot compute FGS countrate & magnitude for a guide star' in str(e_info.value)


@pytest.mark.skip(reason="Removed test until this functionality is added back into the GSSS and this code")
def test_bad_sdss_gz_limits():
    """Test that when SDSSgMag and SDSSzMag color differences are bad, that conversion is skipped"""

    gs_id = 'N13I000018'
    guider = 1
    fgs = FGSCountrate(guide_star_id=gs_id, guider=guider)

    # Reset data to a set of constant, fake data with SDSS_g and z having bad color ranges and tmass missing
    fgs.gsc_series = BASE_DATA
    fgs.gsc_series['tmassJMag'] = -999
    fgs.gsc_series['tmassHMag'] = -999
    fgs.gsc_series['tmassKsMag'] = -999
    fgs.gsc_series['SDSSgMag'] = 20
    fgs.gsc_series['SDSSzMag'] = 14

    # Convert to JHK magnitudes
    fgs.j_mag, fgs.j_mag_err, fgs.h_mag, fgs.h_mag_err, fgs.k_mag, fgs.k_mag_err = \
        fgs.calc_jhk_mag(fgs.gsc_series)

    # Check that the conversion method is the next SDSS one that doesn't include SDSS_g-z
    assert fgs.j_convert_method == 'convert_sdssrz_to_jhk'
    assert fgs.h_convert_method == 'convert_sdssrz_to_jhk'
    assert fgs.k_convert_method == 'convert_sdssrz_to_jhk'


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
    assert np.isclose(j_ser, input_j, 1e-5)
    assert np.isclose(h_ser, input_h, 1e-5)
    assert np.isclose(k_ser, input_k, 1e-5)

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

    # Check tuple and series input produces same result
    assert j_tup == j_ser
    assert h_tup == h_ser
    assert k_tup == k_ser

    # Check conversion function matches hand-check here
    assert np.isclose(j_ser, 11.191999999999998, 1e-5)
    assert np.isclose(h_ser, 11.041, 1e-5)
    assert np.isclose(k_ser, 10.749, 1e-5)

    # Check uncertainties
    assert np.isclose(j_err_ser, 0.27265120293170503, 1e-5)
    assert np.isclose(h_err_ser, 0.24635741034709235, 1e-5)
    assert np.isclose(k_err_ser, 0.2439674603454265, 1e-5)


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

    # Check tuple and series input produces same result
    assert j_tup == j_ser
    assert h_tup == h_ser
    assert k_tup == k_ser

    # Check conversion function matches hand-check here
    assert np.isclose(j_ser, 13.029000000000002, 1e-5)
    assert np.isclose(h_ser, 12.331249999999999, 1e-5)
    assert np.isclose(k_ser, 12.613999999999999, 1e-5)

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

    # Check tuple and series input produces same result
    assert j_tup == j_ser
    assert h_tup == h_ser
    assert k_tup == k_ser

    # Check conversion function matches hand-check here
    assert np.isclose(j_ser, 19.419, 1e-5)
    assert np.isclose(h_ser, 32.059, 1e-5)
    assert np.isclose(k_ser, 24.261999999999997, 1e-5)

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

    # Check tuple and series input produces same result
    assert j_tup == j_ser
    assert h_tup == h_ser
    assert k_tup == k_ser

    # Check conversion function matches hand-check here
    assert np.isclose(j_ser, 11.15, 1e-5)
    assert np.isclose(h_ser, 11.67, 1e-5)
    assert np.isclose(k_ser, 11.73, 1e-5)

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

    # Check tuple and series input produces same result
    assert j_tup == j_ser
    assert h_tup == h_ser
    assert k_tup == k_ser

    # Check conversion function matches hand-check here
    assert np.isclose(j_ser, 11.13, 1e-5)
    assert np.isclose(h_ser, 11.75, 1e-5)
    assert np.isclose(k_ser, 11.899999999999999, 1e-4)

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

    # Check tuple and series input produces same result
    assert j_tup == j_ser
    assert h_tup == h_ser
    assert k_tup == k_ser

    # Check conversion function matches hand-check here
    assert np.isclose(j_ser, 10.02, 1e-5)
    assert np.isclose(h_ser, 11.01, 1e-5)
    assert np.isclose(k_ser, 10.99, 1e-5)

    # Check uncertainties
    assert np.isclose(j_err_ser, 0.3268272426848776, 1e-5)
    assert np.isclose(h_err_ser, 0.387910756252002, 1e-5)
    assert np.isclose(k_err_ser, 0.44586387576927555, 1e-5)
