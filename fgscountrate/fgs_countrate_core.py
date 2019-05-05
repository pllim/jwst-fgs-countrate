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

from . import conversions
from . import utils


class FGS_Countrate():
    """
    Class to support conversion from inputting a guide star ID and
    returning the expected FGS count rate. The main method in this
    class is get_fgs_countrate() which does that conversion.

    Parameters
    ----------
    guide_star_id : str
        The ID of the guide star of interest. This corresponds to the
        HST ID input in the Guide Star Catalog 2

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
        self._wave_list = None
        self._throughput_list = None
        self._signal = None

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

        planck = 6.625e-27

        # Set effective wavelengths (in microns)
        eff_wave_dict = {
            'tmassJmag': 1.25,
            'tmassHmag': 1.65,
            'tmassKsMag': 2.17,
            'SDSSgMag': 0.468,
            'SDSSiMag': 0.748,
            'SDSSzMag': 0.8932,
            'JpgMag': 0.466,
            'FpgMag': 0.645,
            'NpgMag': 0.85,
        }

        # Convert to AB mag
        ab_mag_dict = {key:utils.convert_to_abmag(value, key) if key in self._present_mags else -999
                       for key, value in self._all_mag_series.iteritems()}

        ab_mag_list = [ab_mag_dict[k] for k, v in OrderedDict(sorted(
                       eff_wave_dict.items(), key=operator.itemgetter(1))).items()
                       if k in ab_mag_dict.keys()]  # sort to match order of increasing wavelength

        # Sort to order of increasing wavelength
        self._wave_list = sorted(eff_wave_dict.values())

        # Convert AB mag to flux (photons/s/m**2/micron)
        flux_func = lambda wave, abmag: 10.0**(-(abmag+48.6)/2.5)*(1.0e4/(planck*wave))  # TODO confirm EQN

        flux_list = [flux_func(wave, ab_mag_list[i]) if ab_mag_list[i] != -999 else -999
                     for i, wave in enumerate(self._wave_list)]

        # Fill in gaps
        for i in range(len(flux_list)):
            if i == 0 and flux_list[i] == -999:  # missing shortest band
                flux = -99  # TODO Decide what to do here
            elif flux_list[i] == -999:  # missing middle band
                flux = flux[i-1] + (flux_list[i+1] - flux_list[i-1]) * \
                                   (self._wave_list[i] - self._wave_list[i-1]) / \
                                   (self._wave_list[i+1] - self._wave_list[i-1])
            else:
                break
            flux_list[i] = flux

        # Add extended wavelengths and compute flux
        extended_waves = [3.0, 4.0, 5.0, 5.5]
        self._wave_list.extend(extended_waves)

        extend_flux_func = lambda last_flux, wave: last_flux * 2.17e0**2 / wave**2
        extended_flux_list = [extend_flux_func(flux_list[-1], w) for w in extended_waves]
        flux_list.extend(extended_flux_list)  # TODO how do i use list comprehension to append an already existing list

        # FGS-Guider + OTE throughput. # TODO Fill in SDSS
        throughput_g1_dict = {
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
        throughput_g2_dict = {
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

        if self.guider == 1:
            thr_dict = throughput_g1_dict
            cr_conversion = 999  # TODO NEED THIS VALUE FOR GUIDER 1
        elif self.guider == 2:
            thr_dict = throughput_g2_dict
            cr_conversion = 1.55
        else:
            raise ValueError("Guider value must be an integer either 1 or 2")

        # Compute number of photons reaching the FGS-Guider detector per second per micron.
        self._throughput_list = list(OrderedDict(sorted(thr_dict.items(), key=operator.itemgetter(0))).values())
        self._signal = np.array(flux_list) * np.array(self._throughput_list) * 25  # flux * throughput * PM area

        # Integrate flux in photons per second per micron over the FGS-Guider wavelength range. Use trapezoid formula.
        trapezoid = np.zeros(9)
        for i in range(len(trapezoid)):
            trapezoid[i] = (self._wave_list[i + 1] - self._wave_list[i]) * (self._signal[i] + self._signal[i + 1]) / 2.0
        electrons = np.sum(trapezoid)

        # Compute FGS-Guider count rate using conversion of # electrons => 1 count
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

        # Computation of the FGS magnitude
        trap1 = np.zeros(9)
        trap2 = np.zeros(9)

        for i in range(len(trap1)):
            trap1[i] = (self._wave_list[i + 1] - self._wave_list[i]) * \
                       (self._signal[i] + self._signal[i + 1]) / 2
            trap2[i] = (self._wave_list[i + 1] - self._wave_list[i]) * \
                       (self._throughput_list[i] + self._throughput_list[i + 1]) / 2

        sum_signal = np.sum(trap1)
        sum_throughput = np.sum(trap2)

        self.fgs_magnitude = -2.5 * np.log10(
            sum_signal / sum_throughput) + 27.98  # TODO: Check this is also true for Guider 1

        return self.fgs_magnitude

