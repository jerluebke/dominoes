#include "../include/domino_chain.hpp"
#include <cmath>
#include <algorithm>

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

/*
 * Params:
 *  &lambdas    :   vector with spacings for which to calculate the velocities
 *  limit       :   memory size for integration (numbers of possible
 *      subintervals)
 *
 * Returns:
 *  double_vec_2d result,
 *   result[0][i]   :   angular velocity for spacing lambdas[i]
 *   result[1][i]   :   transversal velocity for spacing lambdas[i]
 *
 */
double_vec_2d DominoChain::make_velocity_array(double_vec& lambdas,
        size_t limit = 1000) const
{
    double_vec_2d result;
    result.resize(2, double_vec());

    params p;
    // wrap theta_dot in lambda function
    GslQuad integrator(
            [this](double theta, params p)
            { return theta_dot(theta, p.eta, p.angular); },
            limit);
    for ( double lambda : lambdas )
    {
        p.eta = _eta(lambda);
        p.angular = intrinsic_angular(lambda);
        result[0].push_back(p.angular);
        // Transversal velocity V = (λ+h)/∫dθ/dθ'
        result[1].push_back(
                (lambda + m_h) / integrator.integrate(p, 0, _psi(lambda)));
    }

    return result;
}

/*
 * Params:
 *  initial_angular :   the angular velocity with which the toppeling domino
 *      chain was initiated
 *  lambda          :   spacing of this domino chain
 *  limt            :   see above
 *
 * Returns:
 *  double_vec_2d result,
 *   result[0][i]   :   spacial coordinate x(i) = i*(λ+h)
 *   result[0][i]   :   angular velocity of the i-th domino
 *   result[1][i]   :   transversal velocity at the i-th domino
 *
 */
double_vec_2d DominoChain::make_velocity_array(double initial_angular,
        double lambda, size_t limit = 1000) const
{
    double_vec_2d result;
    result.resize(3, double_vec());

}

double DominoChain::intrinsic_angular(double lambda) const
{
    return 1.0;
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


//
// Private Methods
//

/*
 * function to be integrated
 * return reciprocal of θ'
 *
 * to prepare this methode for GslQuad, wrap it in some functor object (or
 * lambda) with the parameters eta and angular being passed as a params struct
 *
 */
double DominoChain::theta_dot(double theta, double eta,
        double angular) const
{
    double k_value = k(theta, eta);
    double P_over_K_value = P_over_K(theta, angular);
    double theta_dot = angular * std::sqrt(
            (k_value / (k_value - 1))
            * (1 - P_over_K_value / k_value) );
    return 1 / theta_dot;
}

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
    return ( lambda + m_h ) / m_L;
}


double DominoChain::P_over_K(double theta, double initial_angular_val) const
{
    return 2 * m_omega * m_omega * ( std::cos(m_phi) - std::cos(theta-m_phi) )
        / ( initial_angular_val * initial_angular_val );
}

double DominoChain::k(double theta_initial, double eta) const
{
    double theta_dot_rel_sum = 0.0;
    double theta_dot_rel_prod;
    double theta_dot_rel_value;
    double theta_i;

    for (int j = 0; j < m_N; ++j)
    {
        theta_i = theta_initial;
        theta_dot_rel_prod = 1.0;
        for (int i = 0; i < j; ++i)
        {
            // θ'_i+1 / θ'_i
            theta_dot_rel_value = theta_dot_rel(theta_i, eta);
            theta_dot_rel_prod *= theta_dot_rel_value * theta_dot_rel_value;
            // θ_i = θ_i+1
            theta_i = std::asin( eta*std::cos(theta_i) - m_h/m_L) + theta_i;
        }
        theta_dot_rel_sum += theta_dot_rel_prod;
    }

    return 1 + theta_dot_rel_sum;
}

double DominoChain::theta_dot_rel(double theta, double eta) const
{
    // θ'_i+1 / θ'_1 = 1 - η * sin(θ_i) / cos(θ_i+1 - θ_i)
    // θ_i+1 = arcsin(η * cos(θ_i) - h/L) + θ_i
    return 1 - eta * std::sin(theta)
        / std::cos( std::asin( eta*std::cos(theta) - m_h/m_L ));
}

