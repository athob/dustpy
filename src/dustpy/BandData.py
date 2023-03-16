#!/usr/bin/env python
"""
Utility to load and use band data from SVO FPS
"""
import warnings
import numpy as np
from astropy import table, units
from astroquery.svo_fps import SvoFps

from . import _PARAMS as _P

__all__ = ['BandData']


class BandData:
    def __init__(self, filterids=None):
        """

        :param filterids:
        """
        self.filterids = filterids

    def __dir__(self):
        return super(BandData, self).__dir__() + self._filter_props_attributes

    def __repr__(self):
        return super(BandData, self).__repr__()

    def __getattr__(self, item):
        if not any(le.isupper() for le in item) and item in self._filter_props_attributes:
            return units.Quantity(
                self._filter_props[self._filter_props.keys()[self._filter_props_attributes.index(item)]])
        else:
            raise AttributeError("'%s' object has no attribute '%s'" % (self.__class__.__name__, item))

    @property
    def filterids(self):
        return self._filterids

    @filterids.setter
    def filterids(self, filterids):
        self._filterids = filterids
        self._load_filter_props()
        self._load_filters()

    @property
    def num_of_bands(self):
        return len(self.filterids)

    @property
    def wavelengths_grid(self):
        return self._stacked_grid_of_key(_P.SVO_LAMBDA_KEY)

    @property
    def transmissions_grid(self):
        return self._stacked_grid_of_key(_P.SVO_TRANSMI_KEY)

    @property
    def transmissions_length(self):
        return np.array(list(map(len, self._transmission_data)))

    @property
    def transmissions_integral(self):
        return self._transmission_integral

    def split_by_band(self, array):
        if isinstance(array, table.Table):
            return array.group_by(np.repeat(np.arange(self.num_of_bands), self.transmissions_length)).groups
        else:
            return np.split(array, np.cumsum(self.transmissions_length)[:-1])

    def integrate_over_band(self, i, array=None):
        if array is None:
            array = self._transmission_data[i][_P.SVO_TRANSMI_KEY]
        return np.trapz(array, self._transmission_data[i][_P.SVO_LAMBDA_KEY])

    def _load_filter_props(self):
        self._facilities, self._instruments = zip(*map(lambda fid: fid.split('/'), self._filterids))
        self._instruments, self._bands = zip(*map(lambda ins: ins.split('.'), self._instruments))
        from_svofps = table.vstack(list(map(SvoFps.get_filter_list, np.unique(self._facilities))))
        from_svofps.add_index(_P.SVO_FILTER_ID)
        self._filter_props = from_svofps.loc[self.filterids]
        self._filter_props_attributes = [key.lower() for key in self._filter_props.keys()]

    def _load_filters(self):
        self._transmission_data = self.num_of_bands * [table.QTable()]
        self._transmission_integral = self.num_of_bands * [0.]
        for i, filterid in enumerate(self.filterids):
            _temp = SvoFps.get_transmission_data(filterid)  # TODO: unit consistency
            if isinstance(_temp[_P.SVO_TRANSMI_KEY].unit, units.core.UnrecognizedUnit):
                warnings.warn(f"UnrecognizedUnit '{_temp[_P.SVO_TRANSMI_KEY].unit.name}' found for filter ID '{filterid}', removing it as a temporary solution")
                _temp[_P.SVO_TRANSMI_KEY].unit = ''
            self._transmission_data[i] = table.QTable(_temp)
            self._transmission_integral[i] = self.integrate_over_band(i)

    def _stacked_grid_of_key(self, key):
        stack = units.Quantity(table.vstack(self._transmission_data)[key])
        return stack.unmasked if hasattr(stack, 'unmasked') else stack


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    filter_ids = ['Spitzer/IRAC.I1', 'Spitzer/IRAC.I2', 'Spitzer/IRAC.I3', 'Spitzer/IRAC.I4', 'Spitzer/MIPS.24mu']
    filter_ids += ['WFIRST/WFI.R062', 'WFIRST/WFI.Z087', 'WFIRST/WFI.Y106', 'WFIRST/WFI.J129']
    band_data = BandData(filterids=filter_ids)
