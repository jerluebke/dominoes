# -*- coding: utf-8 -*-

from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.string cimport string

ctypedef vector[double] double_vec
ctypedef vector[vector[double]] double_vec_2d
ctypedef vector[result] result_vec


cdef extern from "GslQuad.hpp":
    ctypedef struct result:
        double error
        int status
        string errormsg


cdef extern from "DominoChain.hpp":
    ctypedef struct domino:
        double height
        double width

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

        double intrinsic_angular(const double,
                                 const double) nogil

        double intrinsic_transversal(const double,
                                     const double,
                                     const bool) nogil

        result_vec& get_full_output() nogil

        void set_pieces_to_be_considered(const int) nogil

