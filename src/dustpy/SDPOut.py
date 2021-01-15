#!/usr/bin/env python
"""
Class definition for SynDustPy output objects
"""
from astropy import units

from . import _PARAMS as _P

__all__ = ['SDPOut']


class SDPOut:
    def __init__(self, band_data, spectra, spectrum):
        """

        :param band_data:
        :param spectra:
        """
        self._band_data = band_data
        self._spectra = spectra
        self._spectrum = spectrum

    def __repr__(self):
        return super(SDPOut, self).__repr__()

    @property
    def band_data(self):
        return self._band_data

    @property
    def spectra(self):
        return self._spectra

    @property
    def spectrum(self):
        return self._spectrum

    @property
    def flux_by_band(self):
        fluxes = self.band_data.num_of_bands * [0.]
        for i, band in enumerate(self.spectra):
            fluxes[i] = self.band_data.integrate_over_band(i,
                                                           band[_P.FLUX_TOT_LABEL + _P.OBSERVED_TAG] /
                                                           band[_P.LAMBDA_LABEL]) / \
                        self.band_data.transmissions_integral[i]
        return units.Quantity(fluxes)


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    pass
