# A modified, cythonized version of some vertex classes of the SHGO
# library by Stefan Endres.
# Original files within the SHGO/scipy library:
#    shgo/_shgo_lib/_vertex.py
#    shgo/_shgo_lib/_complex.py
# Modifications by Malte J. Ziebarth (mjz.science@fmvkb.de)
#
# Copyright (C) 2017 Stefan Endres,
#               2023 Malte J. Ziebarth
#
# This code is licensed under the MIT license (see LICENSE).
# SPDX-License-Identifier: MIT

import numpy as np
cimport numpy as cnp
from cython.operator cimport dereference as deref
from libcpp.vector cimport vector
from libcpp.map cimport map as stdmap
from libcpp cimport bool
from libcpp.memory cimport shared_ptr, make_shared
from cython.parallel cimport prange
from libcpp.set cimport set as stdset
cimport cython

cdef extern from "vertex.hpp" namespace "shgopatch" nogil:

    cdef cppclass CppVertexKey "shgopatch::ScalarVertex::VertexKey":
        long hash
        vector[double] x
        CppVertexKey()
        bool operator==(const CppVertexKey& other) const

    cdef cppclass ScalarVertex:
        ScalarVertex()
        ScalarVertex(long, vector[double]) except+
        size_t hash() const
        CppVertexKey key() const
        void set_fval(double)
        double get_fval() const
        void connect(const ScalarVertex&) except+
        void disconnect(const ScalarVertex&) except+
        bool minimiser() const
        bool maximiser() const
        vector[CppVertexKey] star() const
        vector[CppVertexKey] nearest_neighbors() const
        const vector[double]& x() const

    cdef cppclass ScalarVertexPtrSet:
        ScalarVertexPtrSet()
        void insert(shared_ptr[ScalarVertex]) except+
        void erase(shared_ptr[ScalarVertex]) except+

    cdef void proc_minimiser_set(ScalarVertexPtrSet&)

    cdef vector[CppVertexKey] merge_neighborhoods(const ScalarVertex&,
                                                  const ScalarVertex&)

cdef dict _global_vertices = dict()


cdef class VertexKey:
    cdef CppVertexKey key

    def __eq__(self, other):
        if not isinstance(other, VertexKey):
            return False
        cdef VertexKey ovk = other
        return self.key == ovk.key

    def __hash__(self):
        return self.key.hash


cdef class VertexScalarField:
    """
    Add homology properties of a scalar field f: R^n --> R associated with
    the geometry built from the VertexBase class
    """

    cdef shared_ptr[ScalarVertex] vtx
    cdef object _hash
    cdef cnp.ndarray x_a
    cdef tuple x
    cdef bool _x_a_set
    cdef dict __dict__

    def __init__(self, tuple x, field=None, nn=None, index=None, field_args=(),
                 g_cons=None, g_cons_args=()):
        """
        Parameters
        ----------
        x : tuple,
            vector of vertex coordinates
        field : callable, optional
            a scalar field f: R^n --> R associated with the geometry
        nn : list, optional
            list of nearest neighbours
        index : int, optional
            index of the vertex
        field_args : tuple, optional
            additional arguments to be passed to field
        g_cons : callable, optional
            constraints on the vertex
        g_cons_args : tuple, optional
            additional arguments to be passed to g_cons

        """
        cdef size_t i,n
        cdef VertexScalarField vsf
        cdef vector[double] xvec
        cdef long h
        self._hash = hash(x)
        h = self._hash
        n = len(x)
        xvec.resize(n)
        for i in range(n):
            xvec[i] = x[i]
        self.vtx = make_shared[ScalarVertex](h, xvec)
        self.x = x

        self.index = index
        self._x_a_set = False
        self.__dict__ = dict()
        self.__dict__['x'] = x
        self.__dict__['index'] = index

        # Connect the nearest neighbors:
        if nn is not None:
            for i in range(len(nn)):
                vsf = nn[i]
                deref(self.vtx).connect(deref(vsf.vtx))

        
        # Store in the global cache:
        cdef VertexKey vk = VertexKey()
        vk.key = deref(self.vtx).key()
        _global_vertices[vk] = self

    def __setattr__(self, str name, object value):
        cdef double f
        if name == "f":
            f = value
            deref(self.vtx).set_fval(f)
        self.__dict__[name] = value

    @cython.boundscheck(False)
    cdef _set_x_a(self):
        cdef double[::1] xa
        cdef size_t n
        if not self._x_a_set:
            n = deref(self.vtx).x().size()
            self._x_a = np.empty(n)
            xa = self._x_a
            with nogil:
                for i in range(n):
                    xa[i] = deref(self.vtx).x()[i]
                self._x_a_set = True

    def __getattribute__(self, str name):
        cdef size_t i
        cdef long h
        cdef vector[CppVertexKey] nn
        cdef VertexKey vk
        if name == "x_a":
            self._set_x_a()
            return self._x_a
        elif name == "nn":
            vk = VertexKey()
            ret = set()
            nn = deref(self.vtx).nearest_neighbors()
            for i in range(nn.size()):
                vk.key = nn[i]
                ret.add(_global_vertices[vk])
            return ret
        try:
            return self.__dict__[name]
        except:
            return object.__getattribute__(self, name)


    def __hash__(self) -> int:
        return self._hash

    def connect(self, VertexScalarField v):
        """Connects self to another vertex object v.

        Parameters
        ----------
        v : VertexBase or VertexScalarField object
        """
        deref(self.vtx).connect(deref(v.vtx))

    def disconnect(self, VertexScalarField v):
        deref(self.vtx).disconnect(deref(v.vtx))

    def minimiser(self) -> bool:
        """Check whether this vertex is strictly less than all its
           neighbours"""
        return deref(self.vtx).minimiser()

    def maximiser(self) -> bool:
        """
        Check whether this vertex is strictly greater than all its
        neighbours.
        """
        return deref(self.vtx).maximiser()
        
    def star(self):
        cdef vector[CppVertexKey] st = deref(self.vtx).star()
        cdef set res = set()
        cdef VertexKey vk = VertexKey()
        cdef size_t i
        for i in range(st.size()):
            vk.key = st[i]
            res.add(_global_vertices[vk])
        return res


cdef class VertexCacheDict:
    cdef dict backend
    cdef ScalarVertexPtrSet vertex_set

    def __cinit__(self):
        self.backend = dict()

    def __getitem__(self, key):
        return self.backend[key]

    def __setitem__(self, key, VertexScalarField val):
        if key in self.backend:
            del self[key]
        self.backend[key] = val
        self.vertex_set.insert(val.vtx)

    def __iter__(self):
        yield from self.backend.__iter__()

    def __delitem__(self, key):
        cdef VertexScalarField val = self.backend[key]
        self.vertex_set.erase(val.vtx)
        del self.backend[key]

    def __len__(self):
        return len(self.backend)

#
# This code is copied from scipy/optimize/_shgo_lib/_vertex.py
# The only change is the line 'self.cache = VertexCacheDict()',
# where the collections.OrderedDict is replaced by the above
# VertexCacheDict.
#
class VertexCacheBase:
    """Base class for a vertex cache for a simplicial complex."""
    def __init__(self):

        self.cache = VertexCacheDict()
        self.nfev = 0  # Feasible points
        self.index = -1

    def __iter__(self):
        for v in self.cache:
            yield self.cache[v]
        return

    def size(self):
        """Returns the size of the vertex cache."""
        return self.index + 1

    def print_out(self):
        headlen = len(f"Vertex cache of size: {len(self.cache)}:")
        print('=' * headlen)
        print(f"Vertex cache of size: {len(self.cache)}:")
        print('=' * headlen)
        for v in self.cache:
            self.cache[v].print_out()

#
# The bottleneck in the VertexCacheField class is the `proc_minimisers`
# method.
#
def proc_minimisers(self):
    """Check for minimisers."""
    cdef VertexCacheDict vcd = self.cache
    proc_minimiser_set(vcd.vertex_set)


#
# The follow up seems to be Complex.vpool.
# This code is adapted from scipy/optimize/_shgo_lib/_complex.py
#
def vpool(self, tuple origin, tuple supremum):
    cdef size_t n = len(origin)
    cdef size_t i,j
    cdef vector[double] vot, vst, bl, bu
    cdef double d
    vot.reserve(n)
    vst.reserve(n)
    for d in origin:
        vot.push_back(d)
    for d in supremum:
        vst.push_back(d)
    # Initiate vertices in case they don't exist
    cdef VertexScalarField vo = self.V[origin]
    cdef VertexScalarField vs = self.V[supremum]

    # Remove origin - supremum disconnect

    # Find the lower/upper bounds of the refinement hyperrectangle
    cdef vector[CppVertexKey] k0, k1
    cdef bool success
    with nogil:
        bl.resize(n)
        bu.resize(n)
        for i in range(n):
            bl[i] = min(vot[i], vst[i])
            bu[i] = max(vot[i], vst[i])

        # Merge the two set of nearest neighbors.
        # The result will be in k0.
        if not vo.vtx or not vs.vtx:
            with gil:
                raise RuntimeError("Nullpointer dereference.")
        k0 = merge_neighborhoods(deref(vo.vtx), deref(vs.vtx))
        
        # Apply the bl and bu conditions.
        for j in range(k0.size()):
            success = True
            for i in range(n):
                if k0[j].x[i] < bl[i] or k0[j].x[i] > bu[i]:
                    success = False
                    break
            if success:
                k1.push_back(k0[j])

    # Obtain the Python vertex classes that correspond to the k1 vertices:
    cdef CppVertexKey cvk
    cdef VertexKey vk = VertexKey()
    vn_pool = set()
    for cvk in k1:
        vk.key = cvk
        vn_pool.add(_global_vertices[vk])

    return vn_pool


#
# TODO: This following code is not significantly faster than the Python code.
#       Probably needs pure C++ methods for accessing the Vertex information
#       from coordinates v1 and v2.
#
#def split_edge(self, tuple v1, tuple v2):
#    cdef VertexKey vk = VertexKey()
#    cdef VertexScalarField vc
#    cdef size_t i,n
#    cdef vector[double] vct_cpp
#    cdef double v1xi
#    n = len(v1)
#
#    # Load and disconnect the two vertices:
#    vk.key.x.resize(n)
#    for i in range(n):
#        vk.key.x[i] = v1[i]
#    vk.key.hash = hash(v1)
#    vc = _global_vertices[vk]
#    cdef shared_ptr[ScalarVertex] vtx1 = vc.vtx
#    for i in range(n):
#        vk.key.x[i] = v2[i]
#    vk.key.hash = hash(v2)
#    vc = _global_vertices[vk]
#    cdef shared_ptr[ScalarVertex] vtx2 = vc.vtx
#    if not vtx1 or not vtx2:
#        raise RuntimeError("Nullptr dereferenced in `split_edge`.")
#
#    with nogil:
#        # Destroy original edge, if it exists:
#        deref(vtx1).disconnect(deref(vtx2))
#
#        # Compute vertex on centre of edge:
#        vct_cpp.resize(n)
#        for i in range(n):
#            v1xi = deref(vtx1).key().x[i]
#            vct_cpp[i] = (deref(vtx2).key().x[i] - v1xi) / 2.0 + v1xi
#
#    vc = self.V[tuple(vct_cpp[i] for i in range(n))]
#    # Connect to original 2 vertices to the new centre vertex
#    deref(vc.vtx).connect(deref(vtx1))
#    deref(vc.vtx).connect(deref(vtx2))
#    return vc