# -*- coding: utf-8 -*-

from libcpp cimport bool
from libcpp.vector cimport vector

ctypedef vector[double] double_vec
ctypedef vector[vector[double]] double_vec_2d


cdef extern from "GslQuad.hpp":
    ctypedef struct result:
        pass


cdef extern from "DominoChain.hpp":
    ctypedef struct domino:
        pass

    cdef cppclass DominoChain:
        DominoChain(const domino&,
                    int,
                    const int,
                    const double,
                    const double) nogil

       double_vec_2d make_velocity_array(const double_vec&,
                                         const double,
                                         const bool) nogil

       double_vec_2d make_velocity_array(const double,
                                         const double,
                                         const int,
                                         const double,
                                         const bool) nogil

       result get_full_output(const int) nogil except +

       void set_pieces_to_be_considered(const int) nogil

