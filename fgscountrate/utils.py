def convert_to_abmag(value, name):
    """
    Convert magnitude to AB magnitude

    Parameters
    ----------
    value : float
        Value of the band
    name : str
        Name of the band as stated in the GSC column name.
        Options are: 2MASS: tmassJmag, tmassHmag, tmassKsMag
        SDSS: SDSSgMag, SDSSiMag, SDSSzMag
        GSC: JpgMag, FpgMag, IpgMag

    """
    # TODO Fill in SDSS
    mag_constants = {
        'tmassJmag': 0.90,
        'tmassHmag': 1.37,
        'tmassKsMag': 1.85,
        'SDSSgMag': 999,
        'SDSSiMag': 999,
        'SDSSzMag': 999,
        'JpgMag': -0.055,
        'FpgMag': 0.24,
        'NpgMag': 0.48,
    }

    abmag = value + mag_constants[name]

    return abmag
