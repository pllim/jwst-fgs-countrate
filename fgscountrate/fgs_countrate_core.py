"""
==================
FGS Countrate Core
==================

A module to find the expected FGS countrate for any JWST guide star

Code by Shannon Osborne. Contact at sosborne@stsci.edu
"""

from collections import OrderedDict
import operator

import numpy as np
import pandas as pd

from . import conversions
from . import utils

PLANCK = 6.625e-27

# FGS-Guider + OTE throughput. # TODO Fill in SDSS
THROUGHPUT_G1 = {
    0.466: 0.0,
    0.468: 999,
    0.645: 0.23,
    0.748: 999,
    0.85: 0.45,
    0.8932: 999,
    1.25: 0.56,
    1.65: 0.64,
    2.17: 0.78,
    3.0: 0.70,
    4.0: 0.75,
    5.0: 0.71,
    5.5: 0.0,
}
THROUGHPUT_G2 = {
    0.466: 0.0,
    0.468: 999,
    0.645: 0.35,
    0.748: 999,
    0.85: 0.64,
    0.8932: 999,
    1.25: 0.74,
    1.65: 0.60,
    2.17: 0.63,
    3.0: 0.43,
    4.0: 0.73,
    5.0: 0.68,
    5.5: 0.04,
}

class FGS_Countrate():
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
        self._present_mags = None

        # Conversion information
        self.j_convert_method, self.h_convert_method, self.k_convert_method = None, None, None
        self.j_mag, self.h_mag, self.k_mag = None, None, None

        self.fgs_countrate = None
        self.fgs_magnitude = None
        self.fgs_countrate_data = None

    def get_fgs_countrate(self):
        """
        Calculate the FGS countrate value for a guide star based on it's
        ID number. This calculation using the following steps
            1) Query the newest guide star catalog with the guide star ID
            2) Use the magnitudes present for that guide star to calculate
                the J, H, and K magnitudes
            3) Use the J, H, and K magnitudes along with the guider number
                to calculate the FGS countrate
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
        self.j_mag, self.h_mag, self.k_mag = self.convert_mag_to_jhk(self.data)

        # Compute FGS countrate
        self.fgs_countrate = self.compute_fgs_countrate()

        return self.fgs_countrate

    def convert_mag_to_jhk(self, data):
        """
        Calculate the J, H, and K magnitudes of thee input data

        Parameters
        ----------
        data : pandas series
            A pd series containing at least the following stellar
            magnitude data:
            GSC2: JpgMag, FpgMag, NpgMag
            2MASS: tmassJmag, tmassHmag,tmassKsMag
            SDSS: SDSSuMag, SDSSgMag, SDSSrMag, SDSSiMag, SDSSzMag

        Returns
        -------
        j, h, k : floats
            Runs another method which returns the J, H, and K
            magnitudes for the star

        """

        # Pull all the magnitudes from the series
        l = ['JpgMag', 'FpgMag', 'NpgMag', 'tmassJmag', 'tmassHmag',
             'tmassKsMag', 'SDSSgMag', 'SDSSiMag', 'SDSSzMag']

        self._all_mag_series = data.loc[l]

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

        j = method_list[0](data=self.data, output_mag='J')
        h = method_list[1](data=self.data, output_mag='H')
        k = method_list[2](data=self.data, output_mag='K')

        return j, h, k

    def compute_fgs_countrate(self):
        """
        Compute the FGS countrate using all available values from the GSC

        Returns
        ------
        fgs_countrate : float
            The FGS countrate for the guide star based on J, H, and K mags
            and the guider number
        """

        # Create initial dataframe
        name = ['tmassJmag', 'tmassHmag', 'tmassKsMag', 'SDSSgMag', 'SDSSiMag',
                'SDSSzMag', 'JpgMag', 'FpgMag', 'NpgMag']
        wave = [1.25, 1.65, 2.17, 0.468, 0.748, 0.8932, 0.466, 0.645, 0.85]
        df = pd.DataFrame(wave, columns=['Wavelength'], index=name)

        # Add magnitudes
        df = pd.concat([df, self._all_mag_series], axis=1, sort=True)
        df = df.rename(columns={0: 'Mag'})

        # Calculate and add ABMagnitudes
        def ab_mag(row):
            if row.name in self._present_mags:
                return utils.convert_to_abmag(row['Mag'], row.name)
            else:
                return -999

        df['ABMag'] = df.apply(lambda row: ab_mag(row), axis=1)

        # Sort to order of increasing wavelength
        df = df.sort_values('Wavelength')

        # Convert AB mag to flux (photons/s/m**2/micron)
        def calc_flux(row):
            if row['ABMag'] != -999:
                return 10.0 ** (-(row['ABMag'] + 48.6) / 2.5) * (1.0e4 / (PLANCK * row['Wavelength']))  # TODO confirm EQN
            else:
                return -999

        df['Flux'] = df.apply(lambda row: calc_flux(row), axis=1)

        # Fill in gaps
        # TODO Not quite sure about this
        for i in range(len(df)):
            if i == 0 and df.at[df.index[i], "Flux"] == -999:  # missing shortest band
                flux = -99  # TODO Decide what to do here
            elif df.at[df.index[i], "Flux"] == -999:  # missing middle band
                flux = df.at[df.index[i - 1], "Flux"] + (df.at[df.index[i + 1], "Flux"] - df.at[
                    df.index[i - 1], "Flux"]) * \
                                                        (df.at[df.index[i], "Wavelength"] - df.at[
                                                            df.index[i - 1], "Wavelength"]) / \
                                                        (df.at[df.index[i + 1], "Wavelength"] - df.at[
                                                            df.index[i - 1], "Wavelength"])
            elif i == len(df) and df.at[df.index[i], "Flux"] == -999:  # missing longest band
                flux = -99  # TODO Decide what to do here
            else:
                break
            df.at[df.index[i], "Flux"] = flux

        # Add extended wavelengths and compute flux
        extended_waves = np.array([3.0, 4.0, 5.0, 5.5])
        mag = np.full(4, np.nan)
        abmag = np.full(4, np.nan)
        flux = df.at[df.index[-1], "Flux"] * 2.17e0 ** 2 / extended_waves ** 2
        in_data = list(map(list, zip(extended_waves, mag, abmag, flux)))
        df2 = pd.DataFrame(in_data, columns=['Wavelength', 'Mag', 'ABMag', 'Flux'],
                           index=['Extend1', 'Extend2', 'Extend3', 'Extend4'])
        df = pd.concat([df, df2], axis=0)

        if self.guider == 1:
            df['Throughput'] = df['Wavelength'].map(THROUGHPUT_G1)
            cr_conversion = 999  # TODO NEED THIS VALUE FOR GUIDER 1
        elif self.guider == 2:
            df['Throughput'] = df['Wavelength'].map(THROUGHPUT_G2)
            cr_conversion = 1.55
        else:
            raise ValueError("Guider value must be an integer either 1 or 2")

        # Compute number of photons reaching the FGS-Guider detector per second per micron.
        def calc_signal(row):
            return row['Flux'] * np.array(row['Throughput']) * 25

        df['Signal'] = df.apply(lambda row: calc_signal(row), axis=1)

        # Integrate flux in photons per second per micron over the FGS-Guider wavelength range. Use trapezoid formula.
        trapezoid = np.zeros(9)
        for i in range(len(trapezoid)):
            trapezoid[i] = (df.at[df.index[i + 1], "Wavelength"] - df.at[df.index[i], "Wavelength"]) \
                           * (df.at[df.index[i], "Signal"] + df.at[df.index[i + 1], "Signal"]) / 2.0
        electrons = np.sum(trapezoid)

        # Compute FGS-Guider count rate using conversion of # electrons => 1 count
        self.fgs_countrate_data = df
        self.fgs_countrate = electrons / cr_conversion

        return self.fgs_countrate

    def compute_fgs_magnitude(self):
        """
        Compute the FGS magnitude using all available values from the GSC

        Returns
        ------
        fgs_magnitude : float
            The FGS magnitude for the guide star based on J, H, and K mags
            and the guider number
            """

        self.compute_fgs_countrate()
        df = self.fgs_countrate_data

        # Computation of the FGS magnitude  # TODO: It just ignores the last value?
        length = len(df) - 1
        trap1 = np.zeros(length)
        trap2 = np.zeros(length)
        for i in range(length):
            trap1[i] = (df.at[df.index[i + 1], "Wavelength"] - df.at[df.index[i], "Wavelength"]) * \
                       (df.at[df.index[i], "Signal"] + df.at[df.index[i + 1], "Signal"]) / 2
            trap2[i] = (df.at[df.index[i + 1], "Wavelength"] - df.at[df.index[i], "Wavelength"]) * \
                       (df.at[df.index[i], "Throughput"] + df.at[df.index[i + 1], "Throughput"]) / 2

        sum_signal = np.sum(trap1)
        sum_throughput = np.sum(trap2)

        self.fgs_magnitude = -2.5 * np.log10(
            sum_signal / sum_throughput) + 27.98  # TODO: Check this is also true for Guider 1

        return self.fgs_magnitude

