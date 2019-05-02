"""
==================
FGS Countrate Core
==================

A module to find the expected FGS countrate for any JWST guide star

Code by Shannon Osborne. Contact at sosborne@stsci.edu
"""

from collections import OrderedDict
import io
import operator
import requests

import numpy as np
import pandas as pd

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

    def query_gsc(self):
        """
        Query the Guide Star Catalog 2.4.1 using the guide star ID.

        Returns
        -------
        data : pandas series
            A pd series containing the line from the GSC 2.4.1
            corresponding this the specific guide star ID

        """

        # Query GSC
        file_format = 'CSV'
        catalog = 'GSC241'
        url = 'http://gsss.stsci.edu/webservices/vo/CatalogSearch.aspx?' \
              'GSC2ID={}&FORMAT={}&CAT={}&'.format(self.id, file_format, catalog)
        request = requests.get(url).content

        # Read in data
        try:
            data_frame = pd.read_csv(io.StringIO(request.decode('utf-8')), skiprows=1)
            data_frame.replace(r'^\s+$', -999, regex=True, inplace=True)
        except pd.errors.EmptyDataError:
            raise NameError("This Guide Star ID does not exist in GSC2.4.1")

        # Check length of data table and turn it from a data frame to a series
        if len(data_frame) == 1:
            self.data = data_frame.iloc[0]
        else:
            # TODO: May do more fixing here later
            raise ValueError("This Guide Star ID points to multiple lines in GSC2.4.1")

        return self.data

    def get_fgs_countrate(self):
        """
        This will be the main method called!
        """

        # Query GSC to get data
        data = self.query_gsc()

        # Convert to JHK magnitudes
        self.j_mag, self.h_mag, self.k_mag = self.convert_mag_to_jhk(data)

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
                (i, '_tmass_to_jhk'),
                ('SDSSgMag, SDSSzMag',  '_sdssgz_to_jhk'),
                ('SDSSgMag, SDSSiMag',  '_sdssgi_to_jhk'),
                ('SDSSiMag, SDSSzMag',  '_sdssiz_to_jhk'),
                ('JpgMag, NpgMag',      '_gsc2bjin_to_jhk'),
                ('FpgMag, NpgMag',      '_gsc2rfin_to_jhk'),
                ('JpgMag, FpgMag',      '_gsc2bjrf_to_jhk'),
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
            method = getattr(self, getattr(self, '{}_convert_method'.format(i[5].lower())), lambda: "Invalid")
            method_list.append(method)

        j = method_list[0]('J')
        h = method_list[1]('H')
        k = method_list[2]('K')

        return j, h, k

    def _tmass_to_jhk(self, mag):
        """
        No conversion needed. 2MASS input is already in J,H,K band

        Parameters
        ----------
        mag : str
            The magnitude you want to convert to. Options
            are 'J', 'H', or 'K'.
        """
        if mag.upper() == 'J':
            j = self.data['tmassJmag']
            return j
        elif mag.upper() == 'H':
            h = self.data['tmassHmag']
            return h
        elif mag.upper() == 'K':
            k = self.data['tmassKsMag']
            return k

    def _sdssgz_to_jhk(self, mag):
        """
        Convert from SDSS_g mag and SDSS_z mag to J,H,K mag

        Parameters
        ----------
        mag : str
            The magnitude you want to convert to. Options
            are 'J', 'H', or 'K'.
        """
        g = self.data['SDSSgMag']
        z = self.data['SDSSzMag']
        if mag.upper() == 'J':
            j = g - 0.59 - 1.54*(g - z) + 0.20*(g - z)**2 - 0.04*(g - z)**3 + 0.002*(g - z)**4
            return j
        elif mag.upper() == 'H':
            h = g - 0.77 - 1.78*(g - z) + 0.08*(g - z)**2 - 0.04*(g - z)**3 + 0.009*(g - z)**4
            return h
        elif mag.upper() == 'K':
            k = g - 0.87 - 1.70*(g - z) + 0.01*(g - z)**2 - 0.07*(g - z)**3 + 0.001*(g - z)**4
            return k

    def _sdssgi_to_jhk(self, mag):
        """
        Convert from SDSS_g mag and SDSS_i mag to J,H,K mag

        Parameters
        ----------
        mag : str
            The magnitude you want to convert to. Options
            are 'J', 'H', or 'K'.
        """
        g = self.data['SDSSgMag']
        i = self.data['SDSSiMag']
        if mag.upper() == 'J':
            j = g - 0.411 - 2.260*(g - i) + 0.826*(g - i)**2 - 0.317*(g - i)**3 + 0.037*(g - i)**4
            return j
        elif mag.upper() == 'H':
            h = g - 0.597 - 2.400*(g - i) + 0.450*(g - i)**2 - 0.078*(g - i)**3 + 0.00025*(g - i)**4
            return h
        elif mag.upper() == 'K':
            k = g - 0.637 - 2.519*(g - i) + 0.568*(g - i)**2 - 0.151*(g - i)**3 + 0.013*(g - i)**4
            return k

    def _sdssiz_to_jhk(self, mag):
        """
        Convert from SDSS_i mag and SDSS_z mag to J,H,K mag

        Parameters
        ----------
        mag : str
            The magnitude you want to convert to. Options
            are 'J', 'H', or 'K'.
        """
        i = self.data['SDSSiMag']
        z = self.data['SDSSzMag']
        if mag.upper() == 'J':
            j = i - 0.794 - 2.839*(i - z) + 3.071*(i - z)**2 - 3.139*(i - z)**3 + 1.164*(i - z)**4
            return j
        elif mag.upper() == 'H':
            h = i - 1.051 - 5.361*(i - z) + 8.398*(i - z)**2 - 7.240*(i - z)**3 + 2.111*(i - z)**4
            return h
        elif mag.upper() == 'K':
            k = i - 1.127 - 5.379*(i - z) + 6.454*(i - z)**2 - 3.499*(i - z)**3 + 0.057*(i - z)**4
            return k

    def _gsc2bjin_to_jhk(self, mag):
        """
        Convert from GSC2_B_J mag and GSC2_I_N mag to J,H,K mag

        Parameters
        ----------
        mag : str
            The magnitude you want to convert to. Options
            are 'J', 'H', or 'K'.
        """
        b_j = self.data['JpgMag']
        i_n = self.data['NpgMag']
        if mag.upper() == 'J':
            j = b_j - 1.30*(b_j - i_n) - 0.15
            return j
        elif mag.upper() == 'H':
            h = b_j + 0.06*(b_j - i_n)**2 - 1.71*(b_j - i_n) - 0.10
            return h
        elif mag.upper() == 'K':
            k = b_j + 0.06*(b_j - i_n)**2 - 1.78*(b_j - i_n) - 0.11
            return k

    def _gsc2rfin_to_jhk(self, mag):
        """
        Convert from GSC2_R_F mag and GSC2_I_N mag to J,H,K mag

        Parameters
        ----------
        mag : str
            The magnitude you want to convert to. Options
            are 'J', 'H', or 'K'.
        """
        r_f = self.data['FpgMag']
        i_n = self.data['NpgMag']
        if mag.upper() == 'J':
            j = r_f + 0.01*(r_f - i_n)**2 - 1.56*(r_f - i_n) - 0.44
            return j
        elif mag.upper() == 'H':
            h = r_f + 0.25*(r_f - i_n)**2 - 2.17*(r_f - i_n) - 0.67
            return h
        elif mag.upper() == 'K':
            k = r_f + 0.28*(r_f - i_n)**2 - 2.35*(r_f - i_n) - 0.73
            return k

    def _gsc2bjrf_to_jhk(self, mag):
        """
        Convert from GSC2_B_J mag and GSC2_R_F mag to J,H,K mag

        Parameters
        ----------
        mag : str
            The magnitude you want to convert to. Options
            are 'J', 'H', or 'K'.
        """
        b_j = self.data['JpgMag']
        r_f = self.data['FpgMag']
        if mag.upper() == 'J':
            j = b_j - 0.39*(b_j - r_f)**2 - 0.96*(b_j - r_f) - 0.55
            return j
        elif mag.upper() == 'H':
            h = b_j - 0.24*(b_j - r_f)**2 - 1.66*(b_j - r_f) - 0.41
            return h
        elif mag.upper() == 'K':
            k = b_j - 0.26*(b_j - r_f)**2 - 1.70*(b_j - r_f) - 0.45
            return k

    def compute_fgs_countrate(self):
        """Compute the FGS countrate using all available values from the GSC"""

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
        fgs_countrate = electrons / cr_conversion

        return fgs_countrate

    def compute_fgs_magnitude(self):
        """Compute the FGS magnitude using all available values from the GSC"""

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

