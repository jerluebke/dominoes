#pragma once

#include "integrator.hpp"
#include <vector>

typedef std::vector<double> vec_double;
typedef struct domino
{
    double height;
    double width;
} domino;


class DominoChain
{
    public:
        DominoChain (domino d, int N, int D, double mu);

        template<typename F>
        static GslQuad<F> make_integrator();

        vec_double make_intrinsic_speed_vec(vec_double lambdas) const;
        vec_double make_speed_at_x_vec(double lambda) const;
        
        int intrinsic_angular(double lambda, double* result);
        int intrinsic_transversal(double lambda, 
                double intrinsic_angular_value, double* result);

        int angular_at_x(int i, double initial_val, double* result);
        int transversal_at_x(int i, double initial_angular_val,
                double* result);

        // function to be integrated
        double theta_dot(double theta, double lambda,
                double initial_angular_val) const;


    private:
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
        double k(double theta, double lambda) const;
        double theta_rel(double theta, double lambda) const;
        double theta_next(double theta, double lambda) const;

};
