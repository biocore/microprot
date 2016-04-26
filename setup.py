#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, microprot development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import re
import ast
from setuptools import find_packages, setup


# version parsing from __init__ pulled from Flask's setup.py
# https://github.com/mitsuhiko/flask/blob/master/setup.py
_version_re = re.compile(r'__version__\s+=\s+(.*)')

with open('microprot/__init__.py', 'rb') as f:
    hit = _version_re.search(f.read().decode('utf-8')).group(1)
    version = str(ast.literal_eval(hit))

classifiers = [
    'Development Status :: 2 - Pre-Alpha',
    'License :: OSI Approved :: BSD License',
    'Environment :: Console',
    'Topic :: Software Development :: Libraries :: Application Frameworks',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Programming Language :: Python',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.3',
    'Programming Language :: Python :: 3.4',
    'Operating System :: Unix',
    'Operating System :: POSIX',
    'Operating System :: MacOS :: MacOS X',
    'Operating System :: Microsoft :: Windows']


description = 'MICROPROT: annotation pipeline for microbial (meta)proteomes'
with open('README.md') as f:
    long_description = f.read()

keywords = 'protein annotation structure',

setup(name='microprot',
      version=version,
      license='BSD',
      description=description,
      long_description=long_description,
      keywords=keywords,
      classifiers=classifiers,
      author="microprot development team",
      author_email="tkosciolek@ucsd.edu",
      maintainer="microprot development team",
      maintainer_email="tkosciolek@ucsd.edu",
      url='http://microbio.me/microprot',
      test_suite='nose.collector',
      packages=find_packages(),
      install_requires=[
          'click >= 6',
          'scikit-bio >= 0.4.0',
          'burrito >= 0.9'
      ],
      extras_require={'test': ["nose", "pep8", "flake8"],
                      'coverage': ["coverage"]})
