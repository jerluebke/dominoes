#include "../include/DominoChain.hpp"
#include <cmath>


#ifndef DEBUG
#define DEBUG 0
#endif

#ifdef DEBUG
#define DPRINT(x) do { std::cerr << x << '\n'; } while (0)
#else
#define DPRINT(x)
#endif


const double G = 9.80665;


//////////////////
// Constructor  //
//////////////////

DominoChain::DominoChain(
        const domino &d,
        int N,
        int limit,
        double epsabs,
        double epsrel )

    : theta_dot_wrapper(
            [this]( const double theta, const params& p ) -> double
            {
                return theta_dot( theta, p.index, p.eta, p.angular );
            })
    , m_limit( limit ), m_epsabs( epsabs ), m_epsrel( epsrel )
    , m_N( N ), m_L( d.height ), m_h( d.width )
    , m_phi( _phi(m_L, m_h) ), m_omega( _omega(m_L, m_phi) )
{
    // TODO: why doesn't this work when GslQuad members are `const`?
    // m_integrator = GslQuad<double_func> ( theta_dot, limit );
    // m_integrator = std::make_unique<GslQuad<double_func>>(theta_dot, limit);
}


//////////////////////
// Public Methods   //
//////////////////////

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
double_vec_2d DominoChain::make_velocity_array(
        const double_vec& lambdas,
        const double mu,
        const bool full_output )
{
    if ( full_output )
        _full_output_vec.resize( 0 );

    double_vec_2d result;
    result.resize( 2, double_vec() );

    params p;
    p.index = m_N;
    GslQuad< double_func > integrator( theta_dot_wrapper, m_limit );

    for ( double lambda : lambdas )
    {
        p.eta = _eta( lambda );
        p.angular = intrinsic_angular(
                p.eta, _theta_hat(lambda), _R(lambda, mu) );
        result[0].push_back( p.angular );
        result[1].push_back(
                intrinsic_transversal( integrator, p, lambda, _psi(lambda),
                    full_output ));
    }

    return result;
}


/*
 * Params:
 *  initial_angular :   the angular velocity with which the toppling domino
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
double_vec_2d DominoChain::make_velocity_array(
        const double initial_angular,
        const double lambda,
        const int number_of_pieces,
        const double mu,
        const bool full_output )
{
    if ( full_output )
        _full_output_vec.resize( 0 );

    double_vec_2d result;
    result.resize(3, double_vec());

    const double psi = _psi(lambda);
    const double theta_hat = _theta_hat(lambda);
    const double R = _R(lambda, mu);
    params p;
    p.index = 0;
    p.eta = _eta(lambda);
    p.angular = initial_angular;
    GslQuad< double_func > integrator( theta_dot_wrapper, m_limit );

    while ( p.index <= number_of_pieces )
    {
        result[0].push_back( p.index * (lambda + m_h) );
        result[1].push_back( p.angular );
        result[2].push_back(
                intrinsic_transversal( integrator, p, lambda, psi,
                    full_output ));
        p.angular = angular_next( p.index, p.angular, p.eta, theta_hat, R );
        ++p.index;
    }

    return result;
}


double DominoChain::intrinsic_angular(
        const double eta, 
        const double theta_hat,
        const double R) const
{
    double k_value = k( m_N, 0.0, eta );
    double return_val = m_omega * std::sqrt( ((k_value-1)/k_value)
            * ( 2 * ( std::cos(m_phi) - std::cos(theta_hat - m_phi)) )
            / (k_value * R * R - k_value + 1) );
	return return_val;
}


/*
 * Transversal velocity V = (λ+h)/∫dθ/dθ'
 */
double DominoChain::intrinsic_transversal(
        GslQuad<double_func>& integrator,
        const params& p,
        const double lambda,
        const double psi,
        const bool full_output )
{
    double result = ( lambda + m_h ) / integrator.integrate( p, 0, psi, 
            m_epsabs, m_epsrel, full_output );

    if ( full_output )
        _full_output_vec.push_back( integrator.get_result() );

    return result;
}


double DominoChain::angular_next(
        int index, 
        double initial_val, 
        double eta,
        double theta_hat, 
        double R ) const
{
    double k_value = k( index, 0.0, eta );
    double P_over_K_value = P_over_K( theta_hat, initial_val );
    return initial_val * std::sqrt( ((k_value-1)/k_value)
            * (1 - P_over_K_value/k_value) ) / R;
}


//////////////////////
// Private Methods  //
//////////////////////

/*
 * function to be integrated
 * return reciprocal of θ'
 *
 * to prepare this methode for GslQuad, wrap it in some functor object (or
 * lambda) with the parameters eta and angular being passed as a params struct
 *
 */
double DominoChain::theta_dot( const double theta, const int index,
        const double eta, const double angular ) const
{
    double k_value = k( index, theta, eta );
    double P_over_K_value = P_over_K( theta, angular );
    double theta_dot = angular * std::sqrt(
            (k_value / (k_value - 1))
            * (1 - P_over_K_value / k_value) );
    // DPRINT( theta << '\t' << 1/theta_dot );
    return 1 / theta_dot;
}


// Static Members

double DominoChain::_phi( const double L, const double h )
{
    return std::atan2( h, L );
}

double DominoChain::_omega( const double L, const double phi )
{
    return std::sqrt( 3 * G * std::cos(phi) / (2 * L) );
}


// Non-static Members

double DominoChain::_psi( const double lambda ) const
{
    return std::asin( lambda/m_L );
}

double DominoChain::_xi( const double psi ) const
{
    return m_L*std::cos( psi );
}

double DominoChain::_R( const double lambda, const double mu ) const
{
    double xi = _xi( _psi(lambda) );
    return 1 + ( xi + mu * lambda ) / ( xi - mu * m_h );
}

double DominoChain::_theta_hat( const double lambda ) const
{
    return std::acos( m_h / (m_h + lambda) );
}

double DominoChain::_eta( const double lambda ) const
{
    return ( lambda + m_h ) / m_L;
}


double DominoChain::P_over_K( const double theta,
        const double initial_angular_val ) const
{
    return -2 * m_omega * m_omega * ( std::cos(m_phi) - std::cos(theta-m_phi) )
        / ( initial_angular_val * initial_angular_val );
}

double DominoChain::k( const int piece_index, const double theta_initial,
        const double eta ) const
{
    double theta_dot_rel_sum = 0.0;
    double theta_dot_rel_prod;
    double theta_dot_rel_value;
    double theta_i;
    int pieces_to_consider = ( piece_index < m_N ) ? piece_index+1 : m_N;
    // DPRINT("Pieces to consider: " << pieces_to_consider);

    for ( int j = 0; j < pieces_to_consider; ++j )
    {
        theta_i = theta_initial;
        theta_dot_rel_prod = 1.0;
        for ( int i = 0; i <= j; ++i )
        {
            // θ'_i+1 / θ'_i
            theta_dot_rel_value = theta_dot_rel( theta_i, eta );
            theta_dot_rel_prod *= theta_dot_rel_value * theta_dot_rel_value;
            // θ_i = θ_i+1
            theta_i = std::asin( eta*std::cos(theta_i) - m_h/m_L ) + theta_i;
        }
        theta_dot_rel_sum += theta_dot_rel_prod;
    }

    // DPRINT("k = " << (1+theta_dot_rel_sum));
    return 1 + theta_dot_rel_sum;
}

double DominoChain::theta_dot_rel( const double theta,
        const double eta ) const
{
    // θ'_i+1 / θ'_1 = 1 - η * sin(θ_i) / cos(θ_i+1 - θ_i)
    // θ_i+1 = arcsin(η * cos(θ_i) - h/L) + θ_i
    return 1 - eta * std::sin( theta )
        / std::cos( std::asin( eta*std::cos(theta) - m_h/m_L ) );
}

