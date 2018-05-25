# -*- coding: utf-8 -*-
"""
Cython Wrapper for c++ class `DominoChain` as part of this years SOWAS project
at Ruhr-Uni Bochum

This module contains:
    class `domino_tuple` (named tuple class from pythons `collections` module)
    class `PyDominoChain`

Usage instructions see `PyDominoChain` docstring

---

author: Jeremiah Lübke
email: jeremiah.luebke@ruhr-uni-bochum.de
date: May 2018
"""

from cppDominoChain cimport *

from collections import namedtuple
from numbers import Number
import numpy as np
cimport numpy as np


domino_tuple = namedtuple( "Domino", "height width" )


cdef class PyDominoChain:
    """
    PyDominoChain(domino_tuple, N=10, limit=100,
                  epsabs=1.49e-8, epsrel=1.49e-8)

    A PyDominoChain-Instance holds various methods describing the dynamics of a
    chain of dominoes.

    Parameters
    ----------
    domino_tuple : an arbitrary object representing a domino piece, which
        needs provides a `height` and a `width` attribute.
        This module holds a `domino_tuple` class as a tempalte for this task
    N : int, number of pieces to be considered when doing the calculations.
        This attribute can also be set with `self.set_pieces_to_be_considered`
    limit : size of integration workspace in terms of double precision
        intervals _(can be ignored)_
    epsabs, epsrel : absolute and relative error tolerance for integration
        algorithm _(can be ignored)_

    Methods
    -------
        (description see in docstring for each method)
    intrinsic_velocities
    velocities_by_position
    intrinsic_angular
    intrinsic_transversal
    make_video
    set_pieces_to_be_considered

    Attributes
    ----------
    result_dict : dict holding additional information about the integration
        process

    Examples
    --------
        (see attached jupyter notebook)

    """
    # this-pointer
    cdef DominoChain* cpp_dc

    _result_dict = { "err"      : np.zeros(0, np.float64),
                     "status"   : np.zeros(0, np.float64),
                     "msg"      : [] }

    def __cinit__(self,
                  d_tuple,
                  int N = 10,
                  int limit = 100,
                  double epsabs = 1.49e-8,
                  double epsrel = 1.49e-8):
        """
        Initialize self. Full signature see class docstring
        """
        cdef domino d_struct
        try:
            d_struct.height = d_tuple.height
            d_struct.width = d_tuple.width
        except AttributeError:
            raise AttributeError(
                "The Domino Object passed as parameter needs to provide a `height` and a `width` attribute.")
        self.cpp_dc = new DominoChain( d_struct, N, limit, epsabs, epsrel )


    def __dealloc__(self):
        """
        free memory of pointer to cpp-instance
        """
        del self.cpp_dc


    cpdef np.ndarray intrinsic_velocities(self,
                                          np.ndarray lambdas,
                                          double mu,
                                          bool full_output = False):
        """
        intrinsic_velocities(lambdas, mu, full_output=False)

        Parameters
        ----------
        lambdas : 1-dim array with different spacings for which to calculate
            the intrinsic velocities for
        mu : coefficient of friction
        full_output : whether to retreive additional information from
            integrator

        Returns
        -------
        2-dim array with shape (2, len(lambdas)) holding the intrinsic angular
            and transversal velocities
        """
        cdef double_vec_2d cpp_result = self.cpp_dc.make_velocity_array(
            lambdas, mu, full_output )

        if full_output:
            self._set_result_dict( self.cpp_dc.get_full_output() )

        return np.array( cpp_result, dtype=np.float64 )


    cpdef np.ndarray velocities_by_position(self,
                                            double initial_angular,
                                            double spacing,
                                            int number_of_pieces,
                                            double mu,
                                            bool full_output = False):
        """
        velocities_by_position(initial_angular, spacing, number_of_pieces,
            mu, full_output=False)

        Parameters
        ----------
        initial_angular : angular veloclity with which the toppeling of the
            chain was initiated
        spacing : distance between to pieces
        number_of_pieces : number of pieces in the domino chain
        mu : coefficient of friction
        full_output : whether to retreive additional information from
            integrator

        Returns
        -------
        2-dim array with shape (3, number_of_pieces) holding the x coordinate,
            the angular and the transversal velocity by position
        """
        cdef double_vec_2d cpp_result = self.cpp_dc.make_velocity_array(
            initial_angular, spacing, number_of_pieces, mu, full_output )

        if full_output:
            self._set_result_dict( self.cpp_dc.get_full_output() )

        return np.array( cpp_result, dtype=np.float64 )


    cpdef np.ndarray velocities_variable_spacing(self,
                                                 double initial_angular,
                                                 np.ndarray np_lambdas,
                                                 double mu,
                                                 bool full_output = False):
        """
        velocities_variable_spacing(initial_angular, lambdas, mu,
            full_output=False)

        Parameters
        ----------
        initial_angular : double, angular veloclity with which the toppeling of
            the chain was initiated
        lambdas : 1-dim array with different spacings for which to calculate
            the intrinsic velocities for
        mu : double, coefficient of friction
        full_output : bool, whether to retreive additional information from
            integrator

        Returns
        -------
        2-dim array with shape (3, number_of_pieces) holding the x coordinate,
            the angular and the transversal velocity by position
        """
        cdef double_vec_2d cpp_result = self.cpp_dc.make_velocity_array(
            initial_angular, np_lambdas, mu, full_output )

        if full_output:
            self._set_result_dict( self.cpp_dc.get_full_output() )

        return np.array( cpp_result, dtype=np.float64 )


    cpdef double intrinsic_angular(self,
                                   double spacing,
                                   double mu):
        """
        intrinsic_angular(spacing, mu)

        Parameters
        ----------
        spacing : distance bewteen two pieces
        mu : coefficient of friction

        Returns
        -------
        float64, intrinsic angular velocity for given configuration
        """
        return self.cpp_dc.intrinsic_angular( spacing, mu )


    cpdef double intrinsic_transversal(self,
                                       double spacing,
                                       double mu,
                                       bool full_output = False,
                                       bool times_only = False):
        """
        intrinsic_transversal(spacing, angular, full_output)

        Parameters
        ----------
        spacing : distance between two pieces
        mu : coefficient of friction
        full_output : whether to retreive additional information from
            integrator

        Returns
        -------
        float64, intrinsic transversal velocity for given configuration
        """
        cdef double result = self.cpp_dc.intrinsic_transversal(spacing,
                                                               mu,
                                                               full_output,
                                                               times_only)
        if full_output:
            self._set_result_dict( self.cpp_dc.get_full_output() )

        return result


    cpdef int make_video(self,
                         str filename_str,
                         spacing,
                         double initial_angular,
                         double mu,
                         int number_of_pieces = 128,
                         double fps = 30,
                         int length = 512,
                         int width = 64):
        """
        make_video(filename, spacing, initial_angular, mu, number_of_pieces,
            fps, length, width)

        Parameters
        ----------
        filename : name of the file to save the video in
        spacing : scalar or numpy array holding the distances between the
            dominoes
        initial_angular : angular velocity with which the toppeling of the
            dominoes was initiated
        mu : coefficient of friction
        number_of_pieces : number of dominoes in the chain
            This needs to be evenly divisable by `length` since the number of
            pixels to represent each piece on the screen is determined herewith.
            If `spacing` is an array, its length is used instead and this
            argument is ignored.
        fps : frames per seconds of video
        length, width: measurement of video in pixel

        Returns
        -------
        int, status code (0: success, -1: failure)
        """
        filename = filename_str.encode("utf-8")
        if isinstance(spacing, Number):
            return self.cpp_dc.make_video(filename,
                                          initial_angular,
                                          spacing,
                                          mu,
                                          number_of_pieces,
                                          fps,
                                          length,
                                          width)
        else:
            try:
                spacing = np.array(spacing, dtype=np.float64)
            except TypeError, ValueError:
                raise ValueError("`spacing` must be a number or a numpy array")

            return self.cpp_dc.make_video(filename,
                                          initial_angular,
                                          spacing,
                                          mu,
                                          fps,
                                          length,
                                          width)


    def set_pieces_to_be_considered(self, int value):
        """
        set_pieces_to_be_considered(value)

        Parameters
        ----------
        value : new number of pieces to be considered while performing the
            calculations

        Returns
        -------
        None
        """
        self.cpp_dc.set_pieces_to_be_considered(value)


    @property
    def result_dict( self ):
        """
        dict holding additional information about the integration process

        Keys
        ----
        "err" : integration errors
        "status" : gsl status codes
        "msg" : gsl error message givin details in case of `nan`-results

        Notes
        -----
        When self is initialized, the values of this dict are holding empty
        arrays. In order to fill it, provide the `full_output` flag when
        calling methods which use integration.
        To retreive these information, type
            self.result_dict
        """
        return self._result_dict


    cdef void _set_result_dict( self, result_vec output_vec ):
        """
        PRIVATE MEMBER

        writes additional information from integration process in
        self.result_dict - IS CALLED INTERNALLY

        Parameters
        ----------
        output_vec : std::vector<result>, where `result` is a c-struct holding
            additional information concerning the integraion process

        Returns
        -------
        None
        """
        cdef int length = output_vec.size()

        for key in ("err", "status"):
            self._result_dict[key].resize( length, refcheck=False )
        self._result_dict["msg"].clear()

        for index in range(length):
            result_struct = output_vec[index]
            self._result_dict["err"][index]     = result_struct.error
            self._result_dict["status"][index]  = result_struct.status
            self._result_dict["msg"].append(result_struct.errormsg)

        print("updated `result_dict` ...\n"
              "type `self.result_dict` to retreive it")
        return

