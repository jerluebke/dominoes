# -*- coding: utf-8 -*-

from cppDominoChain cimport result, domino, DominoChain

from collections import namedtuple
import numpy as np

ctypedef np.ndarray[ double, ndim=1, mode="c" ] 1d_array_f
ctypedef np.ndarray[ double, ndim=2, mode="c" ] 2d_array_f
ctypedef np.ndarray[ str, ndim=1, mode="c" ] 1d_array_s

domino_tuple = namedtuple( "Domino", "height width" )


cdef class PyDominoChain:
    """
    TODO:
        finish and test implementation
        add usage instructions...
    """
    # this-pointer
    cdef DominoChain cpp_dc

    _result_dict = { "err"      : np.zeros(0, double),
                     "status"   : np.zeros(0, double),
                     "msg"      : np.zeros(0, str) }

    cdef __cinit__( self,
                    domino_tuple d_tuple,
                    int N = 10,
                    int limit = 100,
                    double epsabs = 1.49e-8,
                    double epsrel = 1.49e-8 ):
        cdef domino d_struct
        d_struct.height = d_tuple.height
        d_struct.width = d_tuple.width
        self.cpp_dc = DominoChain( d_struct, N, limit, epsabs, epsrel )


    cdef 2d_array_f intrinsic_velocities( self,
                                          1d_array_f np_lambdas not None,
                                          double mu,
                                          bool full_output = False ):
        """
        TODO
        """
        cdef double_vec cpp_lambdas ( &lambdas[0], &lambdas[lambdas.size-1] )
        cdef double_vec cpp_result = self.cpp_dc.make_velocity_array(
            cpp_lambdas, mu, full_output )
        cdef 2d_array_f np_result = np.array( cpp_result, dtype=double )

        if full_output:
            self._set_result_dict( self.cpp_dc.get_full_output() )

        return np_result


    cdef 2d_array_f velocities_by_position( self,
                                            double initial_angular,
                                            double spacing,
                                            int number_of_pieces,
                                            double mu,
                                            bool full_output ):
        """
        TODO
        """
        cdef double_vec cpp_result = self.cpp_dc.make_velocity_array(
            initial_angular, spacing, number_of_pieces, mu, full_output )
        cdef 2d_array_f np_result = np.array ( cpp_result, dtype=double )

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
        cdef int length = full_output_vec.size()

        for elem in self._result_dict.values():
            elem.resize( length, refcheck=False )

        cdef int index = 0
        for result_struct in full_output_vec:
            self._result_dict["err"][index]     = result_struct.error
            self._result_dict["status"][index]  = result_struct.status
            self._result_dict["msg"][index]     = result_struct.errormsg
            index += 1

        print( "updated `result_dict` ...\n"
               "type `<this-instance>.result_dict` to retreive it")
        return

