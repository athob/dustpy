#!/usr/bin/env python
"""
Package parameters
"""
import os
from astropy.units import Unit
from astropy.cosmology import default_cosmology

__all__ = ['DUSTY_REPO', 'V_REF_LAMBDA', 'LAMBDA_GRID_RATIO', 'LAMBDA_LABEL',
           'LAMBDA_UNIT', 'LAMBDA_LIMITS', 'FLUX_TOT_LABEL', 'FLUX_INP_LABEL',
           'FLUX_UNIT', 'OUTPUT_RESULTS_SPLIT', 'SVO_LAMBDA_KEY',
           'SVO_TRANSMI_KEY', 'SVO_FILTER_ID', 'DEFAULT_COSMO',
           'IN_SHELL_LABEL', 'THICKNESS_SPLIT', 'OBSERVED_TAG']

# DustPy parameters
DUSTY_REPO = os.path.join(os.path.split(os.path.abspath(__file__))[0], 'data', 'dusty')
VERBOSE = 0
LAMBDA_GRID_RATIO = 1.05
LAMBDA_LABEL = 'lambda'
LAMBDA_UNIT = Unit('micron')
LAMBDA_LIMITS = [1.0e-2, 3.6e4] * LAMBDA_UNIT  # TODO
V_REF_LAMBDA = 0.55 * LAMBDA_UNIT
FLUX_TOT_LABEL = 'fTot'
FLUX_INP_LABEL = 'fInp'
FLUX_UNIT = Unit('W/m2')
OUTPUT_RESULTS_SPLIT = 'RESULTS:'

# BandData parameters
SVO_LAMBDA_KEY = 'Wavelength'
SVO_TRANSMI_KEY = 'Transmission'
SVO_FILTER_ID = 'filterID'

# SynDustPy parameters
DEFAULT_COSMO = default_cosmology
IN_SHELL_LABEL = 'r1'
THICKNESS_SPLIT = 'relative thickness:'
OBSERVED_TAG = 'Obs'

# SDPOut parameters
