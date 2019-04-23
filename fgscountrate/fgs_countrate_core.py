"""
==================
FGS Countrate Core
==================

A module to find the expected FGS countrate for any JWST guide star

Code by Shannon Osborne. Contact at sosborne@stsci.edu
"""

from collections import OrderedDict
import io
import os
import requests
import sys

import numpy as np
import pandas as pd


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

    """

    def __init__(self, guide_star_id):
        self.id = guide_star_id

        self.data = None
        self.j_mag, self.h_mag, self.k_mag = None, None, None
        self.fgs_countrate = None

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
        TBD
        """

        # Query GSC to get data
        data = self.query_gsc()

        # Convert to JHK magnitudes
        self.j_mag, self.h_mag, self.k_mag = self.convert_mag_to_jhk(data)

        # call convert_jhk_to_countrate
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

        present_mags = data.loc[l]

        # Dictionary of convert methods
        switcher = OrderedDict([
            ('tmassJmag, tmassHmag, tmassKsMag', '_convert_tmass_jhk'),
            ('SDSSgMag, SDSSzMag',  '_convert_sdssgz_jhk'),
            ('SDSSgMag, SDSSimag',  '_convert_sdssgi_jhk'),
            ('SDSSiMag, SDSSzMag',  '_convert_sdssiz_jhk'),
            ('JpgMag, NpgMag',      '_convert_gsc2bjin_jhk'),
            ('FpgMag, NpgMag',      '_convert_gsc2rfin_jhk'),
            ('JpgMag, FpgMag',      '_convert_gsc2bjrf_jhk'),
        ])

        # Pull the first entry in the OrderedDict that matches what values are present.
        for key, value in switcher.items():
            key_list = key.split(', ')
            if set(key_list).issubset(present_mags):
                name = key
                break
            else:
                raise ValueError('There is not enough information on this guide star to get its J, H and K data')

        # Get the method
        method_name = switcher.get(name, "nothing")
        method = getattr(self, method_name, lambda: "Invalid")

        return method()

    def _convert_tmass_jhk(self):
        """No conversion needed. 2MASS input is already in j,h,k bands"""
        j = self.data['tmassJmag']
        h = self.data['tmassHmag']
        k = self.data['tmassKsmag']
        return j, h, k

    def _convert_sdssgz_jhk(self):
        """Convert from SDSS_g mag and SDSS_z mag to J,H,K mag"""
        g = self.data['SDSSgMag']
        z = self.data['SDSSzMag']
        j = g - 0.59 - 1.54*(g - z) + 0.20*(g - z)**2 - 0.04*(g - z)**3 + 0.002*(g - z)**4
        h = g - 0.77 - 1.78*(g - z) + 0.08*(g - z)**2 - 0.04*(g - z)**3 + 0.009*(g - z)**4
        k = g - 0.87 - 1.70*(g - z) + 0.01*(g - z)**2 - 0.07*(g - z)**3 + 0.001*(g - z)**4
        return j, h, k

    def _convert_sdssgi_jhk(self):
        """Convert from SDSS_g mag and SDSS_i mag to J,H,K mag"""
        g = self.data['SDSSgMag']
        i = self.data['SDSSiMag']
        j = g - 0.411 - 2.260*(g - i) + 0.826*(g - i)**2 - 0.317*(g - i)**3 + 0.037*(g - i)**4
        h = g - 0.597 - 2.400*(g - i) + 0.450*(g - i)**2 - 0.078*(g - i)**3 + 0.00025*(g - i)**4
        k = g - 0.637 - 2.519*(g - i) + 0.568*(g - i)**2 - 0.151*(g - i)**3 + 0.013*(g - i)**4
        return j, h, k

    def _convert_sdssiz_jhk(self):
        """Convert from SDSS_i mag and SDSS_z mag to J,H,K mag"""
        i = self.data['SDSSiMag']
        z = self.data['SDSSzMag']
        j = i - 0.794 - 2.839*(i - z) + 3.071*(i - z)**2 - 3.139*(i - z)**3 + 1.164*(i - z)**4
        h = i - 1.051 - 5.361*(i - z) + 8.398*(i - z)**2 - 7.240*(i - z)**3 + 2.111*(i - z)**4
        k = i - 1.127 - 5.379*(i - z) + 6.454*(i - z)**2 - 3.499*(i - z)**3 + 0.057*(i - z)**4
        return j, h, k

    def _convert_gsc2bjin_jhk(self):
        """Convert from GSC2_B_J mag and GSC2_I_N mag to J,H,K mag"""
        b_j = self.data['JpgMag']
        i_n = self.data['NpgMag']
        j = b_j - 1.30*(b_j - i_n) - 0.15
        h = b_j + 0.06*(b_j - i_n)**2 - 1.71*(b_j - i_n) - 0.10
        k = b_j + 0.06*(b_j - i_n)**2 - 1.78*(b_j - i_n) - 0.11
        return j, h, k

    def _convert_gsc2rfin_jhk(self):
        """Convert from GSC2_R_F mag and GSC2_I_N mag to J,H,K mag"""
        r_f = self.data['FpgMag']
        i_n = self.data['NpgMag']
        j = r_f + 0.01*(r_f - i_n)**2 - 1.56*(r_f - i_n) - 0.44
        h = r_f + 0.25*(r_f - i_n)**2 - 2.17*(r_f - i_n) - 0.67
        k = r_f + 0.28*(r_f - i_n)**2 - 2.35*(r_f - i_n) - 0.73
        return j, h, k

    def _convert_gsc2bjrf_jhk(self):
        """Convert from GSC2_B_J mag and GSC2_R_F mag to J,H,K mag"""
        b_j = self.data['JpgMag']
        r_f = self.data['FpgMag']
        j = b_j - 0.39*(b_j - r_f)**2 - 0.96*(b_j - r_f) - 0.55
        h = b_j - 0.24*(b_j - r_f)**2 - 1.66*(b_j - r_f) - 0.41
        k = b_j - 0.26*(b_j - r_f)**2 - 1.70*(b_j - r_f) - 0.45
        return j, h, k

    def convert_jhk_to_countrate(self):
        # ...
        # self.fgs_countrate = value
        return self.fgs_countrate
