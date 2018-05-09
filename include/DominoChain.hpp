#pragma once

#include "GslQuad.hpp"
#include <vector>

typedef std::vector<double> double_vec;
typedef std::vector< std::vector<double> > double_vec_2d;

typedef struct domino_struct
{
    double height;
    double width;
} domino;

typedef struct params_struct
{
    double eta;
    double angular;
} params;


class DominoChain
{
    public:
        DominoChain (domino d, int N, int D,
                int limit = 100,
                double epsabs = 1.49e-8,
                double epsrel = 1.49e-8);

        double_vec_2d make_velocity_array(double_vec& lambdas,
                double mu) const;
        double_vec_2d make_velocity_array(double initial_angular,
                double lambda, double mu) const;

        double intrinsic_angular(double eta, double theta_hat,
                double R) const;
        // double intrinsic_transversal(double lambda,
        //         double intrinsic_angular_value) const;
        double angular_next(double initial_val, double eta,
                double theta_hat, double R) const;

		void set_relevant(int value)
		{ m_N = value; }

		void set_dominoes(int value)
		{ m_D = value; }

        void set_limit(int value)
        { m_limit = value; }

        void set_epsabs(double value)
        { m_epsabs = value; }

        void set_epsrel(double value)
        { m_epsrel = value; }


    private:
        // function to be integrated
        double theta_dot(double theta, double eta, double angular) const;

        // integration parameters
        int m_limit;
        double m_epsabs;
        double m_epsrel;

        int m_N;          // number of dominoes to be considered
        int m_D;          // total number of dominoes
        // const double m_mu;      // coefficient of friction
        const double m_L;       // domino height
        const double m_h;       // domino width
        double m_phi;           // = arctan(h/L)
        double m_omega;         // eigenfrequency

        double _psi(double lambda) const;       // angle of impact
        double _xi(double psi) const;           // height of impact
        double _R(double lambda,
                double mu) const;               // transmission of angular speed
        double _theta_hat(double lambda) const; // final angle
        double _eta(double lambda) const;       // = (λ+h)\L

        double P_over_K(double theta, double intitial_angular_speed) const;
        double k(int piece_index, double theta_initial, double eta) const;
        double theta_dot_rel(double theta, double eta) const;

};
