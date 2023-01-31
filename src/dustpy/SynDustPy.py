#!/usr/bin/env python
"""
Synthetic photometry with DustPy and SVO FPS
"""
from astropy import units

from . import _PARAMS as _P
from .BandData import BandData
from .DustPy import DustPy
from .SDPOut import SDPOut

__all__ = ['SynDustPy']


class SynDustPy:
    def __init__(self, filterids=None, src_redshift=0., distance=None):
        """

        :param filterids:
        :param src_redshift:
        :param distance: in Mpc
        """
        self._band_data = BandData(filterids)
        self.src_redshift = src_redshift
        self.distance = distance

    def __repr__(self):
        return super(SynDustPy, self).__repr__()

    @property
    def band_data(self):
        return self._band_data

    @property
    def src_redshift(self):
        return self._src_redshift

    @src_redshift.setter
    def src_redshift(self, src_redshift):
        self._src_redshift = src_redshift
        self._build_dustpy()

    @property
    def distance(self):
        return self._distance

    @distance.setter
    def distance(self, distance):
        if distance is None:
            if self.src_redshift:
                distance = _P.DEFAULT_COSMO.get().lookback_distance(self.src_redshift).to(units.Mpc)
            else:
                raise ValueError("Either :param src_redshift or :param distance should be greater than 0")
        if not isinstance(distance, units.Quantity):
            distance = units.Quantity(distance, unit=units.Mpc)
        self._distance = distance

    def run(self, **kwargs):
        spectrum, results, inputs, full_spectrum = self._dustpy.run(with_full_spectrum=True, **kwargs)
        expansion = 1+self.src_redshift
        spectrum[_P.LAMBDA_LABEL] *= expansion
        full_spectrum[_P.LAMBDA_LABEL] *= expansion
        inner_shell = results[_P.IN_SHELL_LABEL]
        outer_shell = inner_shell * float(inputs.split(_P.THICKNESS_SPLIT)[1].split('\n')[0])
        from_inner = (inner_shell / self.distance).si ** 2
        from_outer = (outer_shell / self.distance).si ** 2
        spectrum[_P.FLUX_TOT_LABEL] *= self.band_data.transmissions_grid  # .unmasked
        spectrum[_P.FLUX_INP_LABEL] *= self.band_data.transmissions_grid  # .unmasked
        spectrum[_P.FLUX_TOT_LABEL+_P.OBSERVED_TAG] = spectrum[_P.FLUX_TOT_LABEL] * from_outer
        spectrum[_P.FLUX_INP_LABEL+_P.OBSERVED_TAG] = spectrum[_P.FLUX_INP_LABEL] * from_inner
        full_spectrum[_P.FLUX_TOT_LABEL+_P.OBSERVED_TAG] = full_spectrum[_P.FLUX_TOT_LABEL] * from_outer
        full_spectrum[_P.FLUX_INP_LABEL+_P.OBSERVED_TAG] = full_spectrum[_P.FLUX_INP_LABEL] * from_inner
        return SDPOut(self.band_data, self.band_data.split_by_band(spectrum), full_spectrum)

    def _build_dustpy(self):
        self._dustpy = DustPy(lambda_grid=self.band_data.wavelengths_grid/(1+self.src_redshift))


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    filter_ids = ['Spitzer/IRAC.I1', 'Spitzer/IRAC.I2', 'Spitzer/IRAC.I3', 'Spitzer/IRAC.I4', 'Spitzer/MIPS.24mu']
    syndustpy = SynDustPy(filterids=filter_ids, distance=778*units.kpc)
    output = syndustpy.run(star_temperature=10000, v_band_optical_depth=1e-1)
