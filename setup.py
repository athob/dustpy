#!/usr/bin/env python
import os

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


def copy_dir():
    dir_path = 'data'
    base_dir = os.path.join('./src/dustpy', dir_path)
    for (dirpath, dirnames, files) in os.walk(base_dir):  # TODO: in dusty only source/, data/ and userpar.inc
        for f in files:
            yield os.path.join(dirpath.split('/', 3)[-1], f)


# metadata are set in the below file, but use this here to avoid warnings.
__author__ = __copyright__ = __credits__ = __license__ = __version__ = __maintainer__ = __email__ = __status__ = None
exec(open(os.path.split(os.path.abspath(__file__))[0] + "/src/dustpy/__metadata__.py").read())

long_description = ""

setup(name='dustpy',
      version=__version__,
      author=__author__,
      author_email=__email__,
      maintainer=__maintainer__,
      maintainer_email=__email__,
      url="https://github.com/athob/dustpy",
      description="",
      long_description=long_description,
      long_description_content_type="text/markdown",
      classifiers=[
          __status__,
          "Environment :: Console",
          "Intended Audience :: Science/Research",
          __license__,
          "Natural Language :: English",
          "Operating System :: Unix",
          "Programming Language :: Python :: 3",
          "Topic :: Database",
          "Topic :: Scientific/Engineering :: Astronomy",
          "Topic :: Software Development :: Version Control :: Git"
      ],
      python_requires='>=3',
      packages=['dustpy'],
      package_dir={'': 'src'},
      package_data={'dustpy': [f for f in copy_dir()]},
      include_package_data=True,
      install_requires=['numpy', 'pandas', 'astropy>=4.2', 'astroquery']  # >=0.4.2']
      )
