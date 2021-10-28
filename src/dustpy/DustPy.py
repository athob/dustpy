#!/usr/bin/env python
"""
Dusty wrapper
"""
import os
import warnings
from glob import glob
import re
import tempfile
import subprocess
from operator import mul

import numpy as np
import pandas as pd
from astropy.units import Quantity, Unit
from astropy.table import QTable

from . import _PARAMS as _P

__all__ = ['DustPy']


class DustPy:
    def __init__(self, lambda_grid=None):
        """

        :param lambda_grid: in microns
        """
        self._build_tempdir = tempfile.TemporaryDirectory()
        self._generate_release()
        self.lambda_grid = lambda_grid

    def __repr__(self):
        return super(DustPy, self).__repr__()

    @property
    def directory_path(self):
        return self._build_tempdir.name + '/dusty'

    @property
    def data_path(self):
        return self.directory_path + '/data'

    @property
    def lambda_grid_path(self):
        return self.data_path + '/lambda_grid.dat'

    @property
    def dusty_path(self):
        return self.directory_path + '/dusty'

    @property
    def lambda_grid(self):
        return Quantity(self._lambda_grid[self._lambda_grid_indices], unit=_P.LAMBDA_UNIT)

    @lambda_grid.setter
    def lambda_grid(self, lambda_grid):
        self._lambda_grid = lambda_grid
        self._clean_release()
        self._prepare_lambda_grid()
        self._build_release()

    def run(self, **kwargs):
        with_full = kwargs.pop('with_full_spectrum', False)
        keep_output = kwargs.pop('keep_output_files', False)  # TODO
        with self._write_input_file(**kwargs) as input_file:
            _temp = subprocess.run([self.dusty_path, input_file.name, _P.VERBOSE],
                                   cwd=self.directory_path, capture_output=True)
            input_file_base = os.path.splitext(input_file.name)[0]
            results, inputs = self._read_output_file(input_file_base + '.out')
            flux_spectrum = self._read_spectrum_file(input_file_base + '.stb')
        if not keep_output:
            for output_file in glob(input_file_base + '*'):
                os.remove(output_file)
        to_return = (flux_spectrum[self._lambda_grid_indices], results, inputs) + with_full*(flux_spectrum,)
        return to_return

    def _generate_release(self):
        # generate release in temp directory
        subprocess.run([_P.DUSTY_REPO + '/../generate_release_light.sh', self._build_tempdir.name])

    def _prepare_lambda_grid(self):
        if self._lambda_grid is None:
            self._lambda_grid = np.loadtxt(self.lambda_grid_path)
            self._lambda_grid_indices = np.arange(self._lambda_grid.shape[0])
        else:
            if isinstance(self._lambda_grid, Quantity):
                self._lambda_grid = self._lambda_grid.to(_P.LAMBDA_UNIT).value.astype('float')
            if type(self._lambda_grid).__module__ != np.__name__:
                self._lambda_grid = np.array(self._lambda_grid)
            _temp_lg_size = self._lambda_grid.shape[0]
            # add reference lambda and limits to grid if absent
            _temp = np.hstack([_P.V_REF_LAMBDA, _P.LAMBDA_LIMITS]).value
            self._lambda_grid = np.hstack([self._lambda_grid, _temp[~np.in1d(_temp, self._lambda_grid)]])
            # sort wavelengths in increasing order
            _sorter = np.argsort(self._lambda_grid)
            self._lambda_grid = self._lambda_grid[_sorter]
            # interpolate between wavelengths with ratio higher than required
            _log10 = np.log10(self._lambda_grid)
            _ext_inds = np.cumsum(np.hstack([0, np.ceil(np.diff(_log10)/np.log10(_P.LAMBDA_GRID_RATIO)).astype('int')]))
            self._lambda_grid = 10 ** np.interp(np.arange(_ext_inds[-1] + 1), _ext_inds, _log10)
            # determine indices that return original grid
            self._lambda_grid_indices = _ext_inds[np.argsort(_sorter)][:_temp_lg_size]
            # write in lambda grid file
            np.savetxt(self.lambda_grid_path, self._lambda_grid, header=' nL = %d' % self._lambda_grid.shape[0])

    def _build_release(self):
        # run makefile in temp directory
        _temp = subprocess.run('make', cwd=self.directory_path, capture_output=True)

    def _clean_release(self):
        pass

    def _write_input_file(self, **kwargs):
        """
        :param kwargs: star_temperature, v_band_optical_depth, inner_edge_temperature
        :return:
        """
        input_file = tempfile.NamedTemporaryFile(mode='w', prefix='dusty', suffix='.inp', dir=self._build_tempdir.name)
        #  ----------------------------------------------------------------------
        #     Input data for DUSTY
        #  ----------------------------------------------------------------------
        #  This is an input file for radiative transfer code DUSTY, version 4.0.
        #  NOTE: this input file is not compatible with old versions of Dusty
        #  due to the added new input options. Examples of each input option are
        #  given at the end of this file. For a more detailed description please
        #  refer to the Manual.
        #  The input file has a free format, text and empty lines can be entered
        #  arbitrarily. All lines that start with the '*' sign are copied to the
        #  output, and can be used to print out notes and comments. This option
        #  can also be useful when the program fails for some mysterious reason
        #  and you want to compare its output with an exact copy of the input line
        #  as it was read in before processing by DUSTY. The occurrence of relevant
        #  numerical input, which is entered in standard FORTRAN conventions, is
        #  flagged by the equal sign `='. The only restrictions are that all required
        #  input entries must be specified, and in the correct order; the most likely
        #  source of an input error is failure to comply with these requirements.
        #  Recall, also, that FORTRAN requires a carriage return termination of the
        #  file's last line if it contains relevant input. Single entries are always
        #  preceded by the equal sign, `=', and must be padded by blanks on both sides;
        #  the terminating blank can be optionally preceded with a comma. For example:
        #  T = 10,000 K as well as Temperature = 1.E4 degrees and simply T = 10000.00
        #  are all equivalent, legal input entries (note that comma separations of long
        #  numbers are permitted).  Some input is entered as a list, in which case the
        #  first member is preceded by `=' and each subsequent member must be preceded
        #  by a blank (an optional comma can be entered before the blank for additional
        #  separation); for example, Temperatures  = 1E4, 2E4 30,000. Because of the
        #  special role of '=' as a flag for input entry, care must be taken not to
        #  introduce any '=' except when required.  All text following the  '%' sign
        #  is ignored (as in TeX) and this can be used to comment out material that
        #  includes '=' signs. For example, different options for the same physical
        #  property may require a different number of input entries. By commenting out
        #  with '%', all options may be retained in the input file with only the
        #  relevant one switched on.
        # >
        # * ----------------------------------------------------------------------
        # * NOTES:
        # * Sample input file (sphere1.inp)
        # * Spherical dust distribution with
        # * constant density profile
        # * heated by external radiation with a Black Body of 5000K
        # * for composite dust grain 70% silicates 30% carbon
        # * for 3 optical depth (30,40,50)
        # * ----------------------------------------------------------------------
        # I. GEOMETRY %(sphere\sphere_matrix\slab)
        #      geometry = sphere
        # II. PHYSICAL PARAMETERS
        #      1) Central source  %(on\off)
        #                 central = on
        #      1.1) Shape: %(black_body\engelkd_marengo\power_law\file_lambda_f_lambda\file_f_lambda\file_f_nu)
        #                 Spectral shape = black_body
        #                 Number of BB = 1
        #                 Temperature = 5800 K
        #      1.2) Scale: %(flux\Lum_r1\energy_den\dilutn_fac\T1)
        #                 Scale = T1   % Td at the inner boundary
        #                 Td = 1200  K
        #      2) External source  %(on\off)
        #                 external = off
        #      3) Dust Properties
        #      3.1 Chemical composition %(common_grain_composite\common_and_addl_grain\tabulated)
        #                 optical properties index = common_grain_composite
        #                 Abundances for supported grain types:
        #                 Sil-Ow  Sil-Oc  Sil-DL  grf-DL  amC-Hn  SiC-Pg
        #            x =  0.00    0.70    0.00    0.00    0.30    0.00
        #                 SIZE DISTRIBUTION = MRN
        #                 Tsub = 1500.
        #      4) Density Distribution %(powd\expd\rdw\rdwa\usr_suppld)
        #                  density type = POWD
        #                  number of powers = 1
        #                  shell's relative thickness = 1000.
        #                  power = 0.
        #      5) Optical Depth: %(linear\logarithmic\user_supplied)
        #                  grid type = logarithmic % log grid
        #                  lambda0 = .55 micron    % fiducial wavelength
        #                  % minimal optical depth @ fiducial wavelength
        #                  tau(min) = 10.0
        #                  % maximal optical depth @ fiducial wavelength
        #                  tau(max) = 1000.0
        #                  number of models = 1
        # ----------------------------------------------------------------------
        # III. NUMERICS
        #   - accuracy for flux conservation = 0.10
        # ----------------------------------------------------------------------
        # IV. OUTPUT PARAMETERS
        #       The flags governing file production are as follows:
        #       If flag.eq.0 the particular file(s) is not produced. If flag.eq.1
        #       all model results are in corresponding files with extensions 'spp'
        #       (sp.properties), 'stb' (spectra), 'itb' (images and visibilities,
        #       if chosen), 'rtb' (radial profiles) and 'mtb' (messages).  If
        #       flag.eq.2 each model result is in a separate corresponding file,
        #       with visibilities contained in fname.i##. If the images flag.eq.3
        #       the visibilities will be in separate files fname.v## (the flag for
        #       visibilities has to be the same as for images).
        #       Note that choosing imaging output requires additional input data
        #       (please refer to the exmaples below or to the Manual).
        #       FILE DESCRIPTION                               FLAG
        #      ------------------------------------------------------------
        #      - detailed spectra for each model;           fname.s### = 2
        #      - images at specified wavelengths;           fname.i### = 0
        #      - en.density at specified radii;             fname.j### = 0
        #      - radial profiles for each model;            fname.r### = 2
        #      - detailed run-time messages;                fname.m### = 2
        #      -------------------------------------------------------------
        # The end of the input parameters listing
        # **********************************************************************
        input_file.file.write("""
I. GEOMETRY
     geometry = %s
""" % ['sphere', 'sphere_matrix', 'slab'][0])
        input_file.file.write("""
II. PHYSICAL PARAMETERS
     1) Central source
                central = %s
""" % ['on', 'off'][0])
        input_file.file.write("""
     1.1) Shape:
                Spectral shape = %s
""" % ['black_body', 'engelkd_marengo', 'power_law', 'file_lambda_f_lambda', 'file_f_lambda', 'file_f_nu'][0])
        input_file.file.write("""
                Number of BB = 1
                Temperature = %f K
""" % kwargs.get('star_temperature', 5800))
        input_file.file.write("""
     1.2) Scale:
                Scale = %s
""" % ['flux', 'Lum_r1', 'energy_den', 'dilutn_fac', 'T1'][-1])
        input_file.file.write("""
                Td = %f K
""" % kwargs.get('inner_edge_temperature', 1200))
        input_file.file.write("""
     2) External source
                external = %s
""" % ['on', 'off'][1])
        input_file.file.write("""
     3) Dust Properties
     3.1 Chemical composition
                optical properties index = %s
""" % ['common_grain_composite', 'common_and_addl_grain', 'tabulated'][0])
        input_file.file.write("""
                Abundances for supported grain types:
                Sil-Ow  Sil-Oc  Sil-DL  grf-DL  amC-Hn  SiC-Pg
           x =  0.00    0.70    0.00    0.00    0.30    0.00
                SIZE DISTRIBUTION = %s
""" % ['mrn', 'modified_mrn', 'kmh'][0])
        input_file.file.write("""
                Tsub = 1500.
     4) Density Distribution
                 density type = %s
""" % ['powd', 'expd', 'rdw', 'rdwa', 'usr_suppld'][0])
        input_file.file.write("""
                 number of powers = 1
                 shell's relative thickness = 1000.
                 power = 0.
     5) Optical Depth:
                 grid type = %s
""" % ['linear', 'logarithmic', 'user_supplied'][1])
        # lambda0  : fiducial wavelength
        # tau(min) : minimal optical depth @ fiducial wavelength
        # tau(max) : maximal optical depth @ fiducial wavelength
        input_file.file.write("""
                 lambda0 = .55 micron
                 tau(min) = %f
                 tau(max) = %f
                 number of models = 1
""" % (2*(kwargs.get('v_band_optical_depth', 10.),)))
        input_file.file.write("""
  III. NUMERICS
     - accuracy for flux conservation = 0.10
  IV. OUTPUT PARAMETERS
        FILE DESCRIPTION                               FLAG
       ------------------------------------------------------------
       - detailed spectra for each model;           fname.s### = 1
       - images at specified wavelengths;           fname.i### = 0
       - en.density at specified radii;             fname.j### = 0
       - radial profiles for each model;            fname.r### = 0
       - detailed run-time messages;                fname.m### = 0
       -------------------------------------------------------------
  The end of the input parameters listing
  **********************************************************************
""")
        input_file.file.flush()
        return input_file

    @staticmethod
    def _read_output_file(path):
        with open(path, 'r') as output_file:
            output = output_file.read()
        try:
            header, results = output.split(_P.OUTPUT_RESULTS_SPLIT)
        except ValueError:
            print(output, flush = True)
            raise ValueError("Couldn't read output file")
        inputs = re.split(r'.inp +', re.split(r'==+', header)[-1])[1]
        names, results = re.split(r'==+', results)[:2]
        names = names.split('\n')[2].split()[1:]
        units = list(map(lambda name: Unit(re.search(r'\((.*)\)', name).group(1) if name[-1] == ')' else '1'), names))
        units[5] = Unit('arcsec')  # TODO: MEH
        names = list(map(lambda name: name.split('(')[0], names))
        results, warning = (results + _P.OUTPUT_RESULTS_WARN).split(_P.OUTPUT_RESULTS_WARN)[:2]
        if warning:
            warnings.warn("IN DUSTY OUTPUT FILE\n" + _P.OUTPUT_RESULTS_WARN + warning)
        results = list(map(float, results.split()[1:]))
        results = QTable(list(zip(map(mul, results, units))), names=names)
        return results, inputs

    @staticmethod
    def _read_spectrum_file(path):
        flux_spectrum = pd.read_fwf(path, skiprows=4)
        flux_spectrum[_P.FLUX_TOT_LABEL] *= flux_spectrum[_P.FLUX_TOT_LABEL][0]
        flux_spectrum[_P.FLUX_INP_LABEL] *= flux_spectrum[_P.FLUX_INP_LABEL][0]
        flux_spectrum.drop(columns='#', index=0, inplace=True)
        units = dict(zip(flux_spectrum.columns, flux_spectrum.shape[1]*[Unit('1')]))
        units[_P.LAMBDA_LABEL] = _P.LAMBDA_UNIT
        units[_P.FLUX_TOT_LABEL] = units[_P.FLUX_INP_LABEL] = _P.FLUX_UNIT
        flux_spectrum = QTable.from_pandas(flux_spectrum, units=units)
        return flux_spectrum


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    dustpy = DustPy(lambda_grid=10 ** (np.random.rand(50) * 5 - 2))
    spectrum = dustpy.run(star_temperature=10000, v_band_optical_depth=0.1)[0]
