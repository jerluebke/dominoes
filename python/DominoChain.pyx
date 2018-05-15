# -*- coding: utf-8 -*-

from cppDominoChain cimport result, domino, DominoChain

from collections import namedtuple

domino_tuple = namedtuple("Domino", "height width")


cdef class PyDominoChain:
    """
    TODO:
        finish and test implementation
        add usage instructions...
    """
    cdef DominoChain cpp_dc

    cdef __cinit__(self,
                   domino_tuple d_tuple,
                   int N = 10,
                   int limit = 100,
                   double epsabs = 1.49e-8,
                   double epsrel = 1.49e-8):
        domino d_struct
        d_struct.height = d_tuple.height
        d_struct.width = d_tuple.width
        self.cpp_dc = DominoChain(d_struct, N, limit, epsabs, epsrel)

    cdef intrinsic_velocities(self):
        """
        TODO
        """

    cdef velocities_by_position(self):
        """
        TODO
        """

    cdef set_pieces_to_be_considered(self, const int value):
        """
        TODO
        """
        self.cpp_dc.set_pieces_to_be_considered(value)
