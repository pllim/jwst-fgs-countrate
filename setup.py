#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst


from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

from setuptools import setup, find_packages

VERSION = '2.0.0'

INSTALL_REQUIRES = [
    'numpy',
    'pytest',
    'requests',
    'pandas',
    'astropy'
]

setup(name='fgscountrate',
      version=VERSION,
      description='JWST FGS countrate estimation',
      long_description='JWST FGS countrate and magnitude estimation and related functionality.',
      classifiers=[
        # 'Development Status :: 3 - Alpha',
        # 'License :: OSI Approved :: MIT License',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Topic :: Scientific/Engineering :: Astronomy'
      ],
      keywords='jwst fgs',
      url='https://github.com/spacetelescope/jwst-fgs-countrate',
      author='Shannon Osborne',
      packages=find_packages(),
      install_requires=INSTALL_REQUIRES,
      include_package_data=True,
      zip_safe=False)
