# Setup script.
#
# Author: Malte J. Ziebarth (mjz.science@fmvkb.de)
#
# Copyright (C) 2023 Malte J. Ziebarth
#
# This code is licensed under the MIT license (see LICENSE).
# SPDX-License-Identifier: MIT


from setuptools import setup
from mebuex import MesonExtension, build_ext

ext = MesonExtension('shgofast.vertex', builddir='builddir')

setup(ext_modules=[ext], cmdclass={'build_ext' : build_ext})
