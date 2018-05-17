# -*- coding: utf-8 -*-
# distutils: language = c++
"""
TODO: docstring
"""

from cppDominoChain cimport *

from collections import namedtuple
import numpy as np
cimport numpy as np

#  ctypedef np.ndarray[ double, ndim=1, mode="c" ] 1d_array_f
#  ctypedef np.ndarray[ double, ndim=2, mode="c" ] 2d_array_f
#  ctypedef np.ndarray[ str, ndim=1, mode="c" ] 1d_array_s

domino_tuple = namedtuple( "Domino", "height width" )


cdef class PyDominoChain:
    """
    TODO:
        finish and test implementation
        add usage instructions...
    """
    # this-pointer
    cdef DominoChain* cpp_dc

    _result_dict = { "err"      : np.zeros(0, np.float64),
                     "status"   : np.zeros(0, np.float64),
                     "msg"      : np.zeros(0, str) }

    def __init__( self,
                    d_tuple,
                    int N = 10,
                    int limit = 100,
                    double epsabs = 1.49e-8,
                    double epsrel = 1.49e-8 ):
        if not isinstance(d_tuple, domino_tuple):
            pass
        cdef domino d_struct
        d_struct.height = d_tuple.height
        d_struct.width = d_tuple.width
        self.cpp_dc = new DominoChain( d_struct, N, limit, epsabs, epsrel )


    def __dealloc__( self ):
        del self.cpp_dc


    cdef np.ndarray intrinsic_velocities( self,
                                          np.ndarray np_lambdas,
                                          double mu,
                                          bool full_output = False ):
        """
        TODO
        """
        cdef double_vec cpp_lambdas = np_lambdas
        cdef double_vec_2d cpp_result = self.cpp_dc.make_velocity_array(
            cpp_lambdas, mu, full_output )
        cdef np.ndarray np_result = np.array( cpp_result, dtype=np.float64 )

        if full_output:
            self._set_result_dict( self.cpp_dc.get_full_output() )

        return np_result


    cdef np.ndarray velocities_by_position( self,
                                            double initial_angular,
                                            double spacing,
                                            int number_of_pieces,
                                            double mu,
                                            bool full_output ):
        """
        TODO
        """
        cdef double_vec_2d cpp_result = self.cpp_dc.make_velocity_array(
            initial_angular, spacing, number_of_pieces, mu, full_output )
        cdef np.ndarray np_result = np.array ( cpp_result, dtype=np.float64 )

        if full_output:
            self._set_result_dict( self.cpp_dc.get_full_output() )

        return np_result


    cdef double intrinsic_angular( self,
                                   double spacing,
                                   double mu ):
        """
        TODO
        """
        return self.cpp_dc.intrinsic_angular( spacing, mu )


    cdef double intrinsic_transversal( self,
                                       double spacing,
                                       double angular,
                                       bool full_output ):
        """
        TODO
        """
        cdef result = self.cpp_dc.intrinsic_transversal( spacing,
                                                         angular,
                                                         full_output)
        if full_output:
            self._set_result_dict( self.cpp_dc.get_full_output() )

        return result



    cdef void set_pieces_to_be_considered(self, int value):
        """
        TODO
        """
        self.cpp_dc.set_pieces_to_be_considered(value)


    @property
    def result_dict( self ):
        return self._result_dict


    cdef void _set_result_dict( self, result_vec output_vec ):
        cdef int length = output_vec.size()

        for elem in self._result_dict.values():
            elem.resize( length, refcheck=False )

        cdef int index = 0
        for result_struct in output_vec:
            self._result_dict["err"][index]     = result_struct.error
            self._result_dict["status"][index]  = result_struct.status
            self._result_dict["msg"][index]     = result_struct.errormsg
            index += 1

        print( "updated `result_dict` ...\n"
               "type `<this-instance>.result_dict` to retreive it")
        return

