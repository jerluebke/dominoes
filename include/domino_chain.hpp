#pragma once

#include "integrator.hpp"
#include <vector>

typedef std::vector<double> double_vec;
typedef std::vector< std::vector<double> > double_vec_2d;

typedef struct domino
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
        DominoChain (domino d, int N, int D, double mu);

        double_vec_2d make_velocity_array(double_vec& lambdas,
                size_t limit) const;
        double_vec_2d make_velocity_array(double initial_angular,
                double lambda, size_t limit) const;

        double intrinsic_angular(double lambda) const;
        int intrinsic_transversal(double lambda,
                double intrinsic_angular_value, double* result);

        int angular_at_x(int i, double initial_val, double* result);
        int transversal_at_x(int i, double initial_angular_val,
                double* result);


    private:
        // function to be integrated
        double theta_dot(double theta, double eta, double angular) const;

        const int m_N;          // number of dominoes to be considered
        const int m_D;          // total number of dominoes
        const double m_mu;      // coefficient of friction
        const double m_L;       // domino height
        const double m_h;       // domino width
        double m_phi;           // = arctan(h/L)
        double m_omega;         // eigenfrequency

        double _psi(double lambda) const;       // angle of impact
        double _xi(double psi) const;           // height of impact
        double _R(double lambda) const;         // transmission of angular speed
        double _theta_hat(double lambda) const; // final angle
        double _eta(double lambda) const;       // = (λ+h)\L

        double P_over_K(double theta, double intitial_angular_speed) const;
        double k(double theta_initial, double eta) const;
        double theta_dot_rel(double theta, double eta) const;

};
