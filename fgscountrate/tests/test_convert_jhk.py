import copy
import itertools

import pandas as pd

from fgscountrate.fgs_countrate_core import FGS_Countrate


# Create data
# Aiming to match output from: fgs = FGS_Countrate(guide_star_id="N13I000018"); data = fgs.query_gsc()
values = ['N13I000018', 420900912, 273.207, 65.5335, 8.30302e-05, 0.000185965, 14.9447, 14.0877,
          13.7468, 13.339, 12.993, 12.901, 14.6547, 14.1443, 14.1067]
index = ['hstID', 'gsc1ID', 'ra', 'dec', 'raErr', 'decErr', 'JpgMag', 'FpgMag', 'NpgMag',
         'tmassJmag', 'tmassHmag', 'tmassKsMag', 'SDSSgMag', 'SDSSiMag', 'SDSSzMag']

BASE_DATA = pd.Series(values, index=index)


def test_convert_mag_to_jhk():
    """
    Using a different technique for choosing the conversion method and
    comparing that output to the technique used in convert_mag_to_jhk.
    This test is done by finding all possible combinations of the
    magnitudes that could be used for the conversion and checking
    that they have the same output with both techniques.

    """
    # Initialize and then re-set the data. No querying done.
    fgs = FGS_Countrate(guide_star_id="N13I000018")
    fgs.data = BASE_DATA

    full_list = ['JpgMag', 'FpgMag', 'NpgMag', 'tmassJmag', 'tmassHmag',
                 'tmassKsMag', 'SDSSgMag', 'SDSSiMag', 'SDSSzMag']

    for L in range(0, len(full_list) + 1):
        for subset in itertools.combinations(full_list, L):
            if subset == ():
                continue

            # Compute conversion
            delete_list = list(set(full_list) - set(subset))
            data2 = copy.copy(fgs.data)
            for i in delete_list:
                data2[i] = -999
            try:
                fgs.convert_mag_to_jhk(data=data2)
            except ValueError:
                error = True

            # Compare output to here
            method_names = []
            for i in ['J', 'H', 'Ks']:
                if 'tmass{}mag'.format(i) in list(subset):
                    method_name_test = "_tmass_to_jhk"

                elif all('SDSSgMag' and 'SDSSzMag') in list(subset):
                    method_name_test = "_sdssgz_to_jhk"
                elif all('SDSSgMag' and 'SDSSimag') in list(subset):
                    method_name_test = "_sdssgi_to_jhk"
                elif all('SDSSiMag' and 'SDSSzMag') in list(subset):
                    method_name_test = "_sdssiz_to_jhk"

                elif all('JpgMag' and 'NpgMag') in list(subset):
                    method_name_test = "_gsc2bjin_to_jhk"
                elif all('FpgMag' and 'NpgMag') in list(subset):
                    method_name_test = "_gsc2rfin_to_jhk"
                elif all('JpgMag' and 'FpgMag') in list(subset):
                    method_name_test = "_gsc2bjrf_to_jhk"

                else:
                    method_name_test = 'cannot_convert_to_jhk'
                method_names.append(method_name_test)

                # If you could compute the JHK mags, check the same conversion method as used
                if error is False:
                    failure_message = 'For input {} and band {}: the test called {} while the convert_mag_to_jhk ' \
                                      'method called {}'.format(list(subset), i, method_name_test,
                                                                getattr(fgs, '{}_convert_method'.format(i[0].lower())))

                    assert method_name_test == getattr(fgs, '{}_convert_method'.format(i[0].lower())), failure_message

            # If you can't compute the JHK mags, check methods agree it's not possible
            if error is True:
                assert 'cannot_convert_to_jhk' in method_names


def test_tmass_to_jhk():
    """
    Check the _tmass_to_jhk method produces the expected result
    which is no change to the input values
    """
    # Create instance and get data
    fgs = FGS_Countrate(guide_star_id="N13I000018")
    fgs.data = BASE_DATA

    # Change data
    input_j = 10
    input_h = 11
    input_k = 12
    fgs.data['tmassJmag'] = input_j
    fgs.data['tmassHmag'] = input_h
    fgs.data['tmassKsmag'] = input_k

    # Run method
    j = fgs._tmass_to_jhk('J')
    h = fgs._tmass_to_jhk('H')
    k = fgs._tmass_to_jhk('K')

    # Check output - there should be no change
    assert j == input_j
    assert h == input_h
    assert k == input_k


def test_sdssgz_to_jhk():
    """Check the _sdssgz_to_jhk method produces the expected result """

    # Create instance and get data
    fgs = FGS_Countrate(guide_star_id="N13I000018")
    fgs.data = BASE_DATA

    # Change data
    input_g = 10
    input_z = 11
    fgs.data['SDSSgMag'] = input_g
    fgs.data['SDSSzMag'] = input_z

    # Run method
    j = fgs._sdssgz_to_jhk('J')
    h = fgs._sdssgz_to_jhk('H')
    k = fgs._sdssgz_to_jhk('K')

    # Do calculation
    val = input_g - input_z
    output_j = input_g - 0.59 - 1.54*val + 0.20*val**2 - 0.04*val**3 + 0.002*val**4
    output_h = input_g - 0.77 - 1.78*val + 0.08*val**2 - 0.04*val**3 + 0.009*val**4
    output_k = input_g - 0.87 - 1.70*val + 0.01*val**2 - 0.07*val**3 + 0.001*val**4

    # Check output
    assert j == output_j
    assert h == output_h
    assert k == output_k


def test_sdssgi_to_jhk():
    """Check the _sdssgi_to_jhk method produces the expected result """

    # Create instance and get data
    fgs = FGS_Countrate(guide_star_id="N13I000018")
    fgs.data = BASE_DATA

    # Change data
    input_g = 10
    input_i = 11
    fgs.data['SDSSgMag'] = input_g
    fgs.data['SDSSiMag'] = input_i

    # Run method
    j = fgs._sdssgi_to_jhk('J')
    h = fgs._sdssgi_to_jhk('H')
    k = fgs._sdssgi_to_jhk('K')

    # Do calculation
    val = input_g - input_i
    output_j = input_g - 0.411 - 2.260*val + 0.826*val**2 - 0.317*val**3 + 0.037*val**4
    output_h = input_g - 0.597 - 2.400*val + 0.450*val**2 - 0.078*val**3 + 0.00025*val**4
    output_k = input_g - 0.637 - 2.519*val + 0.568*val**2 - 0.151*val**3 + 0.013*val**4

    # Check output
    assert j == output_j
    assert h == output_h
    assert k == output_k


def test_sdssiz_to_jhk():
    """Check the _sdssiz_to_jhk method produces the expected result """

    # Create instance and get data
    fgs = FGS_Countrate(guide_star_id="N13I000018")
    fgs.data = BASE_DATA

    # Change data
    input_i = 10
    input_z = 11
    fgs.data['SDSSiMag'] = input_i
    fgs.data['SDSSzMag'] = input_z

    # Run method
    j = fgs._sdssiz_to_jhk('J')
    h = fgs._sdssiz_to_jhk('H')
    k = fgs._sdssiz_to_jhk('K')

    # Do calculation
    val = input_i - input_z
    output_j = input_i - 0.794 - 2.839*val + 3.071*val**2 - 3.139*val**3 + 1.164*val**4
    output_h = input_i - 1.051 - 5.361*val + 8.398*val**2 - 7.240*val**3 + 2.111*val**4
    output_k = input_i - 1.127 - 5.379*val + 6.454*val**2 - 3.499*val**3 + 0.057*val**4

    # Check output
    assert j == output_j
    assert h == output_h
    assert k == output_k


def test_gsc2bjin_to_jhk():
    """Check the _gsc2bjin_to_jhk method produces the expected result """

    # Create instance and get data
    fgs = FGS_Countrate(guide_star_id="N13I000018")
    fgs.data = BASE_DATA

    # Change data
    input_bj = 10
    input_in = 11
    fgs.data['JpgMag'] = input_bj
    fgs.data['NpgMag'] = input_in

    # Run method
    j = fgs._gsc2bjin_to_jhk('J')
    h = fgs._gsc2bjin_to_jhk('H')
    k = fgs._gsc2bjin_to_jhk('K')

    # Do calculation
    val = input_bj - input_in
    output_j = input_bj - 1.30*val - 0.15
    output_h = input_bj + 0.06*val**2 - 1.71*val - 0.10
    output_k = input_bj + 0.06*val**2 - 1.78*val - 0.11

    # Check output
    assert j == output_j
    assert h == output_h
    assert k == output_k


def test_gsc2rfin_to_jhk():
    """Check the _gsc2rfin_to_jhk method produces the expected result """

    # Create instance and get data
    fgs = FGS_Countrate(guide_star_id="N13I000018")
    fgs.data = BASE_DATA

    # Change data
    input_rf = 10
    input_in = 11
    fgs.data['FpgMag'] = input_rf
    fgs.data['NpgMag'] = input_in

    # Run method
    j = fgs._gsc2rfin_to_jhk('J')
    h = fgs._gsc2rfin_to_jhk('H')
    k = fgs._gsc2rfin_to_jhk('K')

    # Do calculation
    val = input_rf - input_in
    output_j = input_rf + 0.01*val**2 - 1.56*val - 0.44
    output_h = input_rf + 0.25*val**2 - 2.17*val - 0.67
    output_k = input_rf + 0.28*val**2 - 2.35*val - 0.73

    # Check output
    assert j == output_j
    assert h == output_h
    assert k == output_k


def test_gsc2bjrf_to_jhk():
    """Check the _gsc2bjrf_to_jhk method produces the expected result """

    # Create instance and get data
    fgs = FGS_Countrate(guide_star_id="N13I000018")
    fgs.data = BASE_DATA

    # Change data
    input_bj = 10
    input_rf = 11
    fgs.data['JpgMag'] = input_bj
    fgs.data['FpgMag'] = input_rf

    # Run method
    j = fgs._gsc2bjrf_to_jhk('J')
    h = fgs._gsc2bjrf_to_jhk('H')
    k = fgs._gsc2bjrf_to_jhk('K')

    # Do calculation
    val = input_bj - input_rf
    output_j = input_bj - 0.39*val**2 - 0.96*val - 0.55
    output_h = input_bj - 0.24*val**2 - 1.66*val - 0.41
    output_k = input_bj - 0.26*val**2 - 1.70*val - 0.45

    # Check output
    assert j == output_j
    assert h == output_h
    assert k == output_k


