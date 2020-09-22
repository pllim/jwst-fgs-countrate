"""
==================
FGS Countrate Core
==================

A module to find the expected FGS countrate and magnitude for any JWST guide star.
Done through
    Querying the Guide Star Catalog
    Converting GSC band data to J, H, and K magnitudes
    Calculating the FGS countrate and FGS magnitude

Code by Shannon Osborne. Contact at sosborne@stsci.edu
"""

from collections import OrderedDict
import copy

import numpy as np
import pandas as pd

from . import conversions
from . import utils

# Universal Constant
PLANCK = 6.625e-27

# GSC Band Information
GSC_BAND_NAMES = ['tmassJmag', 'tmassHmag', 'tmassKsMag',
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


class FGSCountrate:
    """
    Class to support the conversion from an input guide star ID and
    guider number to the expected FGS count rate, FGS magnitude,
    and their uncertainties.

    Parameters
    ----------
    guide_star_id : str
        The ID of the guide star of interest. This corresponds to the
        HST ID input in the Guide Star Catalog

    guider : int
        The guider number, either 1 or 2

    Example
    -------
    fgs = FGSCountrate(guide_star_id='N13I000018', guider=1)
    cr, cr_err, mag, mag_err = fgs.query_fgs_countrate_magnitude()

    """

    def __init__(self, guide_star_id, guider):
        # Star information
        self.id = guide_star_id
        self.guider = guider
        self.gsc_series = None

        # Queried band information
        self._all_queried_mag_series = None
        self._all_queried_mag_err_series = None
        self._present_queried_mags = None

        # Calculated band information (with JHK overwritten)
        self._all_calculated_mag_series = None
        self._all_calculated_mag_err_series = None
        self._present_calculated_mags = None

        # Conversion information
        self.j_convert_method, self.h_convert_method, self.k_convert_method = None, None, None
        self.j_mag, self.j_mag_err, self.h_mag, self.h_mag_err, self.k_mag, self.k_mag_err, = \
            None, None, None, None, None, None

        self.fgs_countrate, self.fgs_countrate_err = None, None
        self.fgs_magnitude, self.fgs_magnitude_err = None, None
        self.band_dataframe = None

    def query_fgs_countrate_magnitude(self, data_frame=None):
        """
        Calculate the FGS countrate and magnitude values for a guide star
        based on it's ID. This calculation uses the following steps:
            1) Query the newest guide star catalog with the given
                guide star ID
            2) Use the magnitudes for the bands present for that
                guide star to calculate the J, H, and K magnitudes
                if they are not already present
            3) Use all the present bands and the guider number to
                calculate the FGS countrate and magnitude

        Returns
        -------
        fgs_countrate : float
            The FGS countrate for the guide star based on J, H, and K mags
            and the guider number
        fgs_countrate_err : float
            Error of the FGS countrate
        fgs_magnitude : float
            The FGS magnitude for the guide star based on J, H, and K mags
            and the guider number
        fgs_magnitude_err : float
            Error of the FGS magnitude
        """

        # Query GSC to get data on the guide star
        if data_frame is None:
            data_frame = utils.query_gsc(gs_id=self.id, catalog='GSC241')

        # Check length of data table and turn it from a dataframe to a series
        if len(data_frame) == 1:
            self.gsc_series = data_frame.iloc[0]
        else:
            raise ValueError("This Guide Star ID points to multiple lines in GSC2.4.1")

        # Convert to JHK magnitudes
        self.j_mag, self.j_mag_err, self.h_mag, self.h_mag_err, self.k_mag, self.k_mag_err = \
            self.calc_jhk_mag()

        # Compute FGS countrate and magnitude
        self.fgs_countrate, self.fgs_countrate_err, \
            self.fgs_magnitude, self.fgs_magnitude_err = self.calc_fgs_cr_mag_and_err()

        return self.fgs_countrate, self.fgs_countrate_err, self.fgs_magnitude, self.fgs_magnitude_err

    def calc_jhk_mag(self, data=None):
        """
        Calculate the J, H, and K magnitudes of the input data

        Parameters
        ----------
        data : pandas series, optional
            A pd series containing at least the following bands:
            GSC2: JpgMag, FpgMag, NpgMag
            2MASS: tmassJmag, tmassHmag, tmassKsMag
            SDSS: SDSSgMag, SDSSiMag, SDSSzMag
            If not passed, will use self.gsc_series

        Returns
        -------
        j, h, k : floats
            Runs another method which returns the J, H, and K
            magnitudes for the star

        """

        # Pull all the magnitudes from the series
        if data is not None:
            self.gsc_series = data

        self._all_queried_mag_series = self.gsc_series.loc[GSC_BAND_NAMES]

        # Pull magnitude errors for each band, and replace missing errors with 2.5% of the magnitude value
        mag_err_list = [self.gsc_series[ind + 'Err'] if self.gsc_series[ind + 'Err'] != -999
                        else self._all_queried_mag_series[ind] * BAND_ERR for ind in self._all_queried_mag_series.index]
        self._all_queried_mag_err_series = pd.Series(mag_err_list, index=self._all_queried_mag_series.index + 'Err')

        # List of the magnitude names that are not fill values in the series
        self._present_queried_mags = list(self._all_queried_mag_series[self._all_queried_mag_series != -999].index)

        # Dictionary of convert methods
        method_list = []
        for i in ['tmassJmag', 'tmassHmag', 'tmassKsMag']:
            switcher = OrderedDict([
                (i, 'convert_tmass_to_jhk'),
                ('SDSSgMag, SDSSzMag',  'convert_sdssgz_to_jhk'),
                ('SDSSgMag, SDSSiMag',  'convert_sdssgi_to_jhk'),
                ('SDSSiMag, SDSSzMag',  'convert_sdssiz_to_jhk'),
                ('JpgMag, NpgMag',      'convert_gsc2bjin_to_jhk'),
                ('FpgMag, NpgMag',      'convert_gsc2rfin_to_jhk'),
                ('JpgMag, FpgMag',      'convert_gsc2bjrf_to_jhk'),
            ])

            # Pull the first entry in the OrderedDict that matches what values are present.
            for key, value in switcher.items():
                key_list = key.split(', ')

                if set(key_list).issubset(self._present_queried_mags):

                    # Check for faint star limits
                    mags = self._all_queried_mag_series[key_list].values
                    if utils.check_band_below_faint_limits(key_list, mags):
                        continue

                    # Set the conversion method
                    setattr(self, '{}_convert_method'.format(i[5].lower()), value)
                    break

            if getattr(self, '{}_convert_method'.format(i[5].lower())) is None:
                raise ValueError('There is not enough information on this guide star to get its {} magnitude'.format(i))

            # Get the method
            method = getattr(conversions, getattr(self, '{}_convert_method'.format(i[5].lower())), lambda: "Invalid")
            method_list.append(method)

        # Create a new series with the edited data (in case uncertainties were replaced)
        edited_data_series = pd.concat([self._all_queried_mag_series, self._all_queried_mag_err_series])

        # Run conversions
        self.j_mag, self.j_mag_err = method_list[0](data=edited_data_series, output_mag='J')
        self.h_mag, self.h_mag_err = method_list[1](data=edited_data_series, output_mag='H')
        self.k_mag, self.k_mag_err = method_list[2](data=edited_data_series, output_mag='K')

        # Create new attribute with updated series
        self._all_calculated_mag_series = copy.deepcopy(self._all_queried_mag_series)
        self._all_calculated_mag_series.loc[['tmassJmag', 'tmassHmag', 'tmassKsMag']] = \
            self.j_mag, self.h_mag, self.k_mag

        self._all_calculated_mag_err_series = copy.deepcopy(self._all_queried_mag_err_series)
        self._all_calculated_mag_err_series.loc[['tmassJmagErr', 'tmassHmagErr', 'tmassKsMagErr']] = \
            self.j_mag_err, self.h_mag_err, self.k_mag_err

        self._present_calculated_mags = self._present_queried_mags + [a for a in
                                                                      ['tmassJmag', 'tmassHmag', 'tmassKsMag']
                                                                      if a not in self._present_queried_mags]

        return self.j_mag, self.j_mag_err, self.h_mag, self.h_mag_err, self.k_mag, self.k_mag_err

    def _calc_fgs_cr_mag(self, to_compute, band_series, guider_throughput, guider_gain, return_dataframe=False):
        """
        Calculate the FGS countrate and/or magnitude based on
        measured guide star data.

        Parameters
        ----------
        to_compute : str
            What FGS value you want to compute and return.
            Options are 'countrate', 'magnitude', or 'both'.
        band_series : Pandas Series
            A Pandas series where the index contains ['tmassJmag',
            'tmassHmag', 'tmassKsMag', 'SDSSgMag', 'SDSSiMag',
            'SDSSzMag', 'JpgMag', 'FpgMag', 'NpgMag'] and the
            values are the magnitude for each of these bands.
        guider_throughput : dict
            A dictionary of the guider (as set in self.guider)
            throughput taken at the central wavelength for each
            of the different bands
        guider_gain : float
            The gain for the guider (as set in self.guider) in
            electrons per 1 count
        return_dataframe : bool
            True or False for if you would like to return a
            dataframe which contains all the calculated
            values for each of the bands including Wavelength,
            Mag, ABMag, Flux, Throughput, and Signal.

        Returns
        -------
        to_return : list
            A list of the different values based on how the "to_compute"
            and "return_dataframe" parameters are set. The possible values
            in this list are (in order):
                FGS countrate
                FGS magnitude
                Dataframe containing the band information

        """
        to_return = []

        # Create initial dataframe
        df = pd.DataFrame(GSC_BAND_WAVELENGTH, columns=['Wavelength'], index=GSC_BAND_NAMES)

        # Add magnitudes
        df = pd.concat([df, band_series], axis=1, sort=True)
        df.columns = ['Wavelength', 'Mag']
        
        # Sort to order of increasing wavelength
        df = df.sort_values(by=['Wavelength'])

        # Calculate and add ABMagnitudes
        def ab_mag(row):
            if row.name in self._present_calculated_mags:
                return utils.convert_to_abmag(row['Mag'], row.name)
            else:
                return -999

        df['ABMag'] = df.apply(lambda row: ab_mag(row), axis=1)

        # Convert AB mag to flux (photons/s/m**2/micron)
        def calc_flux(row):
            if row['ABMag'] != -999:
                return 10.0 ** (-(row['ABMag'] + 48.6) / 2.5) * (1.0e4 / (PLANCK * row['Wavelength']))
            else:
                return -999

        df['Flux'] = df.apply(lambda row: calc_flux(row), axis=1)

        # Add extended wavelengths and compute flux
        extended_waves = np.array([3.0, 4.0, 5.0, 5.5])
        mag = np.full(4, np.nan)
        abmag = np.full(4, np.nan)
        flux = df.at[df.index[-1], "Flux"] * 2.17**2 / extended_waves**2
        in_data = list(map(list, zip(extended_waves, mag, abmag, flux)))
        df2 = pd.DataFrame(in_data, columns=['Wavelength', 'Mag', 'ABMag', 'Flux'],
                           index=['Extend1', 'Extend2', 'Extend3', 'Extend4'])
        df = pd.concat([df, df2], axis=0)

        # Add throughput data to dataframe
        df['Throughput'] = df['Wavelength'].map(guider_throughput)

        # Compute number of photons reaching the FGS-Guider detector per second per micron.
        def calc_signal(row):
            if row['Flux'] != -999:
                return row['Flux'] * row['Throughput'] * 25
            else:
                return -999

        df['Signal'] = df.apply(lambda row: calc_signal(row), axis=1)

        # Can't compute if we only have J, H, and K  # TODO will be changed later
        if {'tmassJmag', 'tmassHmag', 'tmassKsMag'} == set(self._present_calculated_mags):
            raise ValueError('Cannot compute FGS countrate & magnitude for a guide star ({}) with only 2MASS '
                             'data'.format(self.id))

        # Reset the shortest/2nd shortest bands to 0 if they are missing
        if df.at['JpgMag', 'Signal'] == -999:
            df.at['JpgMag', 'Signal'] = 0
        if df.at['SDSSgMag', 'Signal'] == -999:
            df.at['SDSSgMag', 'Signal'] = 0

        # Throw out any bands for which you don't have data
        df_short = df[df['Signal'] != -999]

        # Integrate flux in photons per second per micron over the FGS-Guider wavelength range. Use trapezoid formula.
        trapezoid = np.zeros(len(df_short) - 1)
        for i in range(len(trapezoid)):
            trapezoid[i] = (df_short.at[df_short.index[i + 1], "Wavelength"] -
                            df_short.at[df_short.index[i], "Wavelength"]) * \
                           (df_short.at[df_short.index[i], "Signal"] +
                            df_short.at[df_short.index[i + 1], "Signal"]) / 2.0
        electrons = np.sum(trapezoid)

        # Compute FGS countrate using conversion of # electrons => 1 count
        fgs_countrate = electrons / guider_gain

        if to_compute.lower() == 'countrate' or to_compute.lower() == 'both':
            to_return.append(fgs_countrate)

        # Compute FGS magnitude
        length = len(df_short) - 1
        trap1 = np.zeros(length)
        trap2 = np.zeros(length)
        for i in range(length):
            trap1[i] = (df_short.at[df_short.index[i + 1], "Wavelength"] -
                        df_short.at[df_short.index[i], "Wavelength"]) * \
                       (df_short.at[df_short.index[i], "Signal"] +
                        df_short.at[df_short.index[i + 1], "Signal"]) / 2.0
            trap2[i] = (df_short.at[df_short.index[i + 1], "Wavelength"] -
                        df_short.at[df_short.index[i], "Wavelength"]) * \
                       (df_short.at[df_short.index[i], "Throughput"] +
                        df_short.at[df_short.index[i + 1], "Throughput"]) / 2.0

        sum_signal = np.sum(trap1)
        sum_throughput = np.sum(trap2)

        if self.guider == 1:
            mag_conversion = MAG_CONVERSION_G1
        elif self.guider == 2:
            mag_conversion = MAG_CONVERSION_G2
        else:
            raise ValueError("Guider value must be an integer either 1 or 2")

        fgs_magnitude = -2.5 * np.log10(sum_signal / sum_throughput) + mag_conversion

        if to_compute.lower() == 'magnitude' or to_compute.lower() == 'both':
            to_return.append(fgs_magnitude)

        if return_dataframe is True:
            to_return.append(df)

        return to_return

    def calc_fgs_cr_mag_and_err(self):
        """
        Compute the FGS countrate and magnitude and their respective errors
        using all available bands from the GSC

        Returns
        -------
        fgs_countrate : float
            The FGS countrate for the guide star based on J, H, and K mags
            and the guider number
        fgs_countrate_err : float
            Error of the FGS countrate
        fgs_magnitude : float
            The FGS magnitude for the guide star based on J, H, and K mags
            and the guider number
        fgs_magnitude_err : float
            Error of the FGS magnitude
        """

        # Set values based on guider
        if self.guider == 1:
            throughput_dict = THROUGHPUT_G1
            cr_conversion = CR_CONVERSION_G1
        elif self.guider == 2:
            throughput_dict = THROUGHPUT_G2
            cr_conversion = CR_CONVERSION_G2
        else:
            raise ValueError("Guider value must be an integer either 1 or 2")

        # Calculate magnitude/countrate
        self.fgs_countrate, self.fgs_magnitude, self.band_dataframe = \
            self._calc_fgs_cr_mag(to_compute='both', band_series=self._all_calculated_mag_series,
                                  guider_throughput=throughput_dict, guider_gain=cr_conversion,
                                  return_dataframe=True)

        # Band Magnitude Error
        cr_err_list = []
        mag_err_list = []
        for band in self._present_calculated_mags:
            band_data_with_err = copy.deepcopy(self._all_calculated_mag_series)
            band_data_with_err[band] += self._all_calculated_mag_err_series[band+'Err']
            cr_band_err, mag_band_err = self._calc_fgs_cr_mag(to_compute='both',
                                                              band_series=band_data_with_err,
                                                              guider_throughput=throughput_dict,
                                                              guider_gain=cr_conversion)
            cr_err_list.append(cr_band_err - self.fgs_countrate)
            mag_err_list.append(mag_band_err - self.fgs_magnitude)

        # Throughput Error - 5%
        new_throughput = {key: val * 1.05 for key, val in throughput_dict.items()}
        cr_tput_err, mag_tput_err = self._calc_fgs_cr_mag(to_compute='both',
                                                          band_series=self._all_calculated_mag_series,
                                                          guider_throughput=new_throughput,
                                                          guider_gain=cr_conversion)
        cr_err_list.append(cr_tput_err - self.fgs_countrate)
        mag_err_list.append(mag_tput_err - self.fgs_magnitude)

        # Gain Error - 5%
        new_gain = cr_conversion * 1.05
        cr_gain_err, mag_gain_err = self._calc_fgs_cr_mag(to_compute='both',
                                                          band_series=self._all_calculated_mag_series,
                                                          guider_throughput=throughput_dict,
                                                          guider_gain=new_gain)
        cr_err_list.append(cr_gain_err - self.fgs_countrate)

        # Integral Error - 5%
        cr_err_list.append(self.fgs_countrate * 0.05)
        mag_err_list.append(self.fgs_magnitude * 0.05)

        # Combine Error
        self.fgs_countrate_err = np.sqrt(sum(i**2 for i in cr_err_list))
        self.fgs_magnitude_err = np.sqrt(sum(i**2 for i in mag_err_list))

        return self.fgs_countrate, self.fgs_countrate_err, self.fgs_magnitude, self.fgs_magnitude_err
