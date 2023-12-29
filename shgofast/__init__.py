# Monkey-patching the Cython code into the SciPy SHGO code.
#
# Copyright (C) 2023 Malte J. Ziebarth
#
# This code is licensed under the MIT license (see LICENSE).
# SPDX-License-Identifier: MIT

from .vertex import VertexScalarField as VSF, proc_minimisers, VertexCacheBase,\
                    vpool #, split_edge
import scipy.optimize._shgo_lib._vertex as shgo_vertex
import scipy.optimize._shgo_lib._complex as shgo_complex

# The function in the Cython code is, for some reason, not bound to instances
# if it is directly inserted into the `VertexCacheField` class.
def _proc_minimisers(self):
    return proc_minimisers(self)

def _vpool(self, origin, supremum):
    return vpool(self, origin, supremum)

#def _split_edge(self, v1, v2):
#    return split_edge(self, v1, v2)


# Save the original class and method:
_original_VertexScalarField = shgo_vertex.VertexScalarField
_original_proc_minimisers = shgo_vertex.VertexCacheField.proc_minimisers
_original_VertexCacheBase = shgo_vertex.VertexCacheBase
_original_vpool = shgo_complex.Complex.vpool
_original_split_edge = shgo_complex.Complex.split_edge

# Monkey patch:
def patch():
    """
    Patches the scipy SHGO routines using 
    """
    shgo_vertex.VertexScalarField = VSF
    shgo_vertex.VertexCacheField.proc_minimisers = _proc_minimisers
    shgo_vertex.VertexCacheBase = VertexCacheBase
    shgo_vertex.VertexCacheField.__bases__ = (VertexCacheBase,)
    shgo_complex.Complex.vpool = _vpool
#    shgo_complex.Complex.split_edge = _split_edge

def unpatch():
    shgo_vertex.VertexScalarField = _original_VertexScalarField
    shgo_vertex.VertexCacheField.proc_minimisers = _original_proc_minimisers
    shgo_vertex.VertexCacheBase = _original_VertexCacheBase
    shgo_vertex.VertexCacheField.__bases__ = (_original_VertexCacheBase,)
    shgo_complex.Complex.vpool = _original_vpool
#    shgo_complex.Complex.split_edge = _original_split_edge
