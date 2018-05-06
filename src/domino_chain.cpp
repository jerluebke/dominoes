#include "../include/domino_chain.hpp"
#include <cmath>

const double G = 9.80665;


//
// Constructor
//

DominoChain::DominoChain(domino d, int N, int D, double mu)
    : m_L(d.height), m_h(d.width), m_N(N), m_D(D), m_mu(mu)
{
    m_phi = std::atan2(m_h, m_L);
    m_omega = std::sqrt(3*G*std::cos(m_phi)/(2*m_L));
}


//
// Public Methods
//

template<typename F>
GslQuad<F> DominoChain::make_integrator()
{

}

vec_double DominoChain::make_intrinsic_speed_vec(vec_double lambdas) const
{

}

vec_double DominoChain::make_speed_at_x_vec(double lambda) const
{

}

int DominoChain::intrinsic_angular(double lambda, double* result)
{

}

int DominoChain::intrinsic_transversal(double lambda,
        double intrinsic_angular_value, double* result)
{

}

int DominoChain::angular_at_x(int i, double initial_val, double* result)
{

}

int DominoChain::transversal_at_x(int i, double initial_angular_val,
        double* result)
{

}

double DominoChain::theta_dot(double theta, double lambda,
        double initial_angular_val) const
{

}


//
// Private Methods
//

double DominoChain::_psi(double lambda) const
{
    return std::asin(lambda/m_L);
}

double DominoChain::_xi(double psi) const
{
    return m_L*std::cos(psi);
}

double DominoChain::_R(double lambda) const
{
    double xi = this->_xi(this->_psi(lambda));
    return 1 + ( xi + m_mu * lambda ) / ( xi - m_mu * m_h );
}

double DominoChain::_theta_hat(double lambda) const
{
    return std::acos( m_h / (m_h + lambda) );
}

double DominoChain::_eta(double lambda) const
{
    return ( lambda + m_h ) / m_h;
}


double DominoChain::P_over_K(double theta, double initial_angular_val) const
{

}

double DominoChain::k(double theta, double lambda) const
{

}

double DominoChain::theta_rel(double theta, double lambda) const
{

}

double DominoChain::theta_next(double theta, double lambda) const
{

}
