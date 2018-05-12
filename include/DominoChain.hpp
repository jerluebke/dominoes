#pragma once

#include "GslQuad.hpp"
#include <vector>
#include <functional>

typedef std::vector<double> double_vec;
typedef std::vector< std::vector<double> > double_vec_2d;

typedef struct domino_struct
{
    double height;
    double width;
} domino;

typedef struct params_struct
{
    int index;
    double eta;
    double angular;
} params;

typedef std::function< double (double, params) > double_func;


class DominoChain
{
    public:
        DominoChain (
                const domino& d,
                int N,
                const int limit = 100,
                const double epsabs = 1.49e-8,
                const double epsrel = 1.49e-8 );

        double_vec_2d make_velocity_array(
                const double_vec& lambdas,
                const double mu ) const;

        double_vec_2d make_velocity_array(
                const double initial_angular,
                const double lambda,
                const int number_of_pieces,
                const double mu ) const;

        double intrinsic_angular(
                const double eta,
                const double theta_hat,
                const double R ) const;

        double intrinsic_transversal(
                const GslQuad<double_func>& integrator,
                const params& p,
                const double lambda,
                const double psi ) const;

        double angular_next(
                const int index,
                const double initial_val,
                const double eta,
                const double theta_hat,
                const double R ) const;

		void set_pieces_to_be_considered( const int value )
		{ m_N = value; }


    private:
        // function to be integrated
        double theta_dot(
                const double theta,
                const int index,
                const double eta,
                const double angular ) const;

        const double_func theta_dot_wrapper;

        // integrator params
        const size_t m_limit;
        const double m_epsabs;
        const double m_epsrel;

        // domino quantities
        int m_N;          // number of dominoes to be considered
        const double m_L;       // domino height
        const double m_h;       // domino width
        const double m_phi;     // = arctan(h/L)
        const double m_omega;   // eigenfrequency

        // used to initialize m_phi and m_omega
        static double _phi( const double L, const double h );
        static double _omega( const double L, const double phi );

        // helper methods
        double _psi( const double lambda ) const;   // angle of impact
        double _xi( const double psi ) const;       // height of impact
        double _R( const double lambda,
                const double mu ) const;            // transmission of angular speed
        double _theta_hat( const double lambda ) const;   // final angle
        double _eta( const double lambda ) const;   // = (λ+h)/L

        double P_over_K(
                const double theta,
                const double intitial_angular_speed ) const;
        double k(
                const int piece_index,
                const double theta_initial,
                const double eta ) const;
        double theta_dot_rel(
                const double theta,
                const double eta ) const;

};
