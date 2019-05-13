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

import numpy as np
import pandas as pd

from . import conversions
from . import utils

# Constants
PLANCK = 6.625e-27

# GSC Band Information
GSC_BAND_NAMES = ['tmassJmag', 'tmassHmag', 'tmassKsMag', 'SDSSgMag', 'SDSSiMag',
                  'SDSSzMag', 'JpgMag', 'FpgMag', 'NpgMag']
GSC_BAND_WAVELENGTH = [1.25, 1.65, 2.17, 0.468, 0.748, 0.8932, 0.466, 0.645, 0.85]

# FGS-Guider + OTE throughput.
THROUGHPUT_G1 = {
    0.466: 0.042,
    0.468: 0.044,
    0.645: 0.487,
    0.748: 0.586,
    0.850: 0.655,
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
    0.466: 0.020,
    0.468: 0.021,
    0.645: 0.390,
    0.748: 0.628,
    0.850: 0.669,
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
MAG_CONVERSION_G1 = -24.7934
MAG_CONVERSION_G2 = -24.7011


class FGS_Countrate:
    """
    Class to support conversion from inputting a guide star ID and
    returning the expected FGS count rate. The main method in this
    class is get_fgs_countrate() which does that conversion.

    Parameters
    ----------
    guide_star_id : str
        The ID of the guide star of interest. This corresponds to the
        HST ID input in the Guide Star Catalog

    guider : int
        The guider number, either 1 or 2

    """

    def __init__(self, guide_star_id, guider):
        # Star information
        self.id = guide_star_id
        self.guider = guider
        self.data = None

        # Band information
        self._all_mag_series = None
        self._all_mag_err_series = None
        self._present_mags = None

        # Conversion information
        self.j_convert_method, self.h_convert_method, self.k_convert_method = None, None, None
        self.j_mag, self.j_mag_err, self.h_mag, self.h_mag_err, self.k_mag, self.k_mag_err, = \
            None, None, None, None, None, None

        self.fgs_countrate, self.fgs_countrate_err = None, None
        self.fgs_magnitude, self.fgs_magnitude_err = None, None
        self.fgs_countrate_data = None  # TODO figure out where to define this

    def get_fgs_countrate_magnitude(self):
        """
        Calculate the FGS countrate and magnitude values for a guide star based on it's
        ID number. This calculation using the following steps
            1) Query the newest guide star catalog with the guide star ID
            2) Use the magnitudes present for that guide star to calculate
                the J, H, and K magnitudes
            3) Use the J, H, and K magnitudes along with the guider number
                to calculate the FGS countrate and magnitude

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

        # Query GSC to get data
        data_frame = utils.query_gsc(gs_id=self.id, catalog='GSC241')

        # Check length of data table and turn it from a data frame to a series
        if len(data_frame) == 1:
            self.data = data_frame.iloc[0]
        else:
            # TODO: May do more fixing here later
            raise ValueError("This Guide Star ID points to multiple lines in GSC2.4.1")

        # Convert to JHK magnitudes
        self.j_mag, self.j_mag_err, self.h_mag, self.h_mag_err, self.k_mag, self.k_mag_err = \
            self.convert_mag_to_jhk(self.data)

        # Compute FGS countrate and magnitude
        self.fgs_countrate, self.fgs_countrate_err, \
        self.fgs_magnitude, self.fgs_magnitude_err = self.compute_fgs_data()

        return self.fgs_countrate, self.fgs_countrate_err, self.fgs_magnitude, self.fgs_magnitude_err

    def convert_mag_to_jhk(self, data):
        """
        Calculate the J, H, and K magnitudes of thee input data

        Parameters
        ----------
        data : pandas series
            A pd series containing at least the following stellar
            magnitude data:
            GSC2: JpgMag, FpgMag, NpgMag
            2MASS: tmassJmag, tmassHmag, tmassKsMag
            SDSS: SDSSgMag, SDSSiMag, SDSSzMag

        Returns
        -------
        j, h, k : floats
            Runs another method which returns the J, H, and K
            magnitudes for the star

        """

        # Pull all the magnitudes from the series
        self._all_mag_series = data.loc[GSC_BAND_NAMES]

        # Pull magnitude errors for each band
        mag_err_list = []
        for ind in self._all_mag_series.index:
            mag_err_list.append(self.data[ind + 'Err'])
        self._all_mag_err_series = pd.Series(mag_err_list, index=self._all_mag_series.index)

        # List of the magnitude names that are not fill values in the series
        self._present_mags = list(self._all_mag_series[self._all_mag_series != -999].index)

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
                if set(key_list).issubset(self._present_mags):
                    setattr(self, '{}_convert_method'.format(i[5].lower()), value)
                    break
            if getattr(self, '{}_convert_method'.format(i[5].lower())) is None:
                raise ValueError('There is not enough information on this '
                                 'guide star to get its {} magnitude'.format(i))

            # Get the method
            method = getattr(conversions, getattr(self, '{}_convert_method'.format(i[5].lower())), lambda: "Invalid")
            method_list.append(method)

        self.j_mag, self.j_mag_err = method_list[0](data=self.data, output_mag='J')
        self.h_mag, self.h_mag_err = method_list[1](data=self.data, output_mag='H')
        self.k_mag, self.k_mag_err = method_list[2](data=self.data, output_mag='K')

        return self.j_mag, self.j_mag_err, self.h_mag, self.h_mag_err, self.k_mag, self.k_mag_err

    def _compute_band_data(self, to_compute, band_data, guider_throughput, guider_gain):

        # Create initial dataframe
        df = pd.DataFrame(GSC_BAND_WAVELENGTH, columns=['Wavelength'], index=GSC_BAND_NAMES)

        # Add magnitudes
        df = pd.concat([df, band_data], axis=1, sort=True)
        df = df.rename(columns={0: 'Mag'})

        # Sort to order of increasing wavelength
        df = df.sort_values(by=['Wavelength'])

        # Calculate and add ABMagnitudes
        def ab_mag(row):
            if row.name in self._present_mags:
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
        flux = df.at[df.index[-1], "Flux"] * 2.17e0 ** 2 / extended_waves ** 2
        in_data = list(map(list, zip(extended_waves, mag, abmag, flux)))
        df2 = pd.DataFrame(in_data, columns=['Wavelength', 'Mag', 'ABMag', 'Flux'],
                           index=['Extend1', 'Extend2', 'Extend3', 'Extend4'])
        df = pd.concat([df, df2], axis=0)

        # Set values based on guider
        df['Throughput'] = df['Wavelength'].map(guider_throughput)
        cr_conversion = guider_gain

        # Compute number of photons reaching the FGS-Guider detector per second per micron.
        def calc_signal(row):
            if row['Flux'] != -999:
                return row['Flux'] * row['Throughput'] * 25
            else:
                return -999

        df['Signal'] = df.apply(lambda row: calc_signal(row), axis=1)

        # Can't compute if we only have J, H, and K  # TODO will be changed
        if set(['tmassJmag', 'tmassHmag', 'tmassKsMag']) == set(self._present_mags):
            raise ValueError('Cannot compute FGS countrate for a guide star ({}) with only 2MASS data'.format(self.id))

        # Reset the shortest/2nd shortest bands to 0 if they are missing
        if df.at['JpgMag', 'Signal'] == -999:
            df.at['JpgMag', 'Signal'] = 0
        if df.at['SDSSgMag', 'Signal'] == -999:
            df.at['SDSSgMag', 'Signal'] = 0

        # # Compute FGS-Guider count rate using conversion of # electrons => 1 count
        # self.fgs_countrate_data = df

        # Throw out any data you don't have
        df_short = df[df['Signal'] != -999]

        # Integrate flux in photons per second per micron over the FGS-Guider wavelength range. Use trapezoid formula.
        trapezoid = np.zeros(len(df_short) - 1)
        for i in range(len(trapezoid)):
            trapezoid[i] = (df_short.at[df_short.index[i + 1], "Wavelength"] -
                            df_short.at[df_short.index[i], "Wavelength"]) * \
                           (df_short.at[df_short.index[i], "Signal"] +
                            df_short.at[df_short.index[i + 1], "Signal"]) / 2.0
        electrons = np.sum(trapezoid)

        fgs_countrate = electrons / cr_conversion

        if to_compute == 'countrate':
            return fgs_countrate

        # Computation of the FGS magnitude
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

        if to_compute == 'magnitude':
            return fgs_magnitude
        elif to_compute == 'both':
            return fgs_countrate, fgs_magnitude

    def compute_fgs_data(self):
        """
        Compute the FGS countrate and magnitude (and their respective errors)
        using all available values from the GSC

        Returns
        ------
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
        self.fgs_countrate, self.fgs_magnitude = self._compute_band_data(to_compute='both',
                                                                         band_data=self._all_mag_series,
                                                                         guider_throughput=throughput_dict,
                                                                         guider_gain=cr_conversion)

        # Calculate uncertainty
        cr_err_list = []
        mag_err_list = []
        for band in self._present_mags:
            band_data_with_err = self._all_mag_series
            band_data_with_err[band] += self._all_mag_err_series[band]
            cr_band_err, mag_band_err = self._compute_band_data(to_compute='both',
                                                                band_data=band_data_with_err,
                                                                guider_throughput=throughput_dict,
                                                                guider_gain=cr_conversion)
            cr_err_list.append(cr_band_err - self.fgs_countrate)
            mag_err_list.append(mag_band_err - self.fgs_magnitude)

        # Throughput Error - 5%
        new_throughput = {key: val * 1.05 for key, val in throughput_dict.items()}
        cr_tput_err, mag_tput_err = self._compute_band_data(to_compute='both',
                                                            band_data=self._all_mag_series,
                                                            guider_throughput=new_throughput,
                                                            guider_gain=cr_conversion)
        cr_err_list.append(cr_tput_err - self.fgs_countrate)
        mag_err_list.append(mag_tput_err - self.fgs_magnitude)

        # Gain Error - 5%
        new_gain = cr_conversion * 1.05
        cr_gain_err, mag_gain_err = self._compute_band_data(to_compute='both',
                                                            band_data=self._all_mag_series,
                                                            guider_throughput=throughput_dict,
                                                            guider_gain=new_gain)
        cr_err_list.append(cr_gain_err - self.fgs_countrate)
        mag_err_list.append(mag_gain_err - self.fgs_magnitude)

        # Integral Error - 5% # TODO - figure out what value to use
        cr_integral_err = self.fgs_countrate * 1.05
        mag_integral_err = self.fgs_magnitude * 1.05
        cr_err_list.append(cr_integral_err - self.fgs_countrate)
        mag_err_list.append(mag_integral_err - self.fgs_magnitude)

        # Combine Error
        self.fgs_countrate_err = np.sqrt(np.sum(i**2 for i in cr_err_list))
        self.fgs_magnitude_err = np.sqrt(np.sum(i**2 for i in mag_err_list))

        return self.fgs_countrate, self.fgs_countrate_err, self.fgs_magnitude, self.fgs_magnitude_err