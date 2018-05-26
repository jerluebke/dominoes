#include "../include/DominoChain.hpp"
#include <cmath>
#include <algorithm>


#ifndef PRINT_EXTRA
#define PRINT_EXTRA 0
#endif

#ifdef PRINT_EXTRA
#define DPRINT(x) do { std::cerr << x << '\n'; } while (0)
#else
#define DPRINT(x)
#endif

#ifndef WITH_R
#define WITH_R 1
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
    , m_phi(_phi( d.height, m_h )), m_omega(_omega( d.height, m_phi ))
{}


//////////////////////
// Public Methods   //
//////////////////////

/*
 * Params:
 *  &lambdas    :   vector with spacings for which to calculate the velocities
 *  mu          :   double, coefficient of friction
 *  full_output :   bool, whether to retreive additional information from
 *      integrator and store it in m_full_output_vec
 *
 * Returns:
 *  double_vec_2d result,
 *   result[0][i]   :   angular velocity for spacing lambdas[i]
 *   result[1][i]   :   transversal velocity for spacing lambdas[i]
 *
 */
const double_vec_2d DominoChain::make_velocity_array(
        const double_vec& lambdas,
        const double mu,
        const bool full_output,
        const bool times_only )
{
    if ( full_output )
        m_full_output_vec.clear();

    double_vec_2d result( times_only ? 1 : 2 );

    params p;
    p.index = m_N;
    GslQuad< double_func > integrator( theta_dot_wrapper, m_limit );

    for ( double lambda : lambdas )
    {
        p.eta = _eta( lambda );
        p.angular = intrinsic_angular(
                p.eta, _theta_hat(lambda), _R(lambda, mu) );
        if ( times_only )
        {
            result[0].push_back(
                    integrator.integrate( p, 0, _psi(lambda),
                        m_epsabs, m_epsrel, full_output ));
        }
        else
        {
            result[0].push_back( p.angular );
            result[1].push_back(
                    transversal( integrator, p, lambda, _psi(lambda),
                        full_output ));
        }
    }

    return result;
}


/*
 * Params:
 *  initial_angular :   double, the angular velocity with which the toppling
 *      domino chain was initiated
 *  lambda          :   double, spacing of this domino chain
 *  mu              :   double, coefficient of friction
 *  number_of_pieces    :   double, number of dominoes in chain
 *  full_output :   bool, whether to retreive additional information from
 *      integrator and store it in m_full_output_vec
 *
 * Returns:
 *  double_vec_2d result,
 *   result[0][i]   :   spacial coordinate x(i) = i*(λ+h)
 *   result[0][i]   :   angular velocity of the i-th domino
 *   result[1][i]   :   transversal velocity at the i-th domino
 *
 */
const double_vec_2d DominoChain::make_velocity_array(
        const double initial_angular,
        const double lambda,
        const int number_of_pieces,
        const double mu,
        const bool full_output,
        const bool times_only )
{
    if ( full_output )
        m_full_output_vec.clear();

    double_vec_2d result( times_only ? 1 : 3 );
    double psi = _psi(lambda);
    const double theta_hat = _theta_hat(lambda);
    const double R = _R(lambda, mu);
    params p;
    p.index = 0;
    p.eta = _eta(lambda);
    p.angular = initial_angular;
    GslQuad< double_func > integrator( theta_dot_wrapper, m_limit );

    while ( p.index <= number_of_pieces )
    {
        if ( p.index == number_of_pieces )
            psi = M_PI/2;

        if ( times_only )
        {
            result[0].push_back(
                    integrator.integrate( p, 0, psi, m_epsabs, m_epsrel,
                        full_output ));
        }
        else
        {
            result[0].push_back( p.index * (lambda + m_h) );
            result[1].push_back( p.angular );
            result[2].push_back(
                    transversal( integrator, p, lambda, psi,
                        full_output ));
        }
        p.angular = angular_next( p.index, p.angular, p.eta, theta_hat, R );
        ++ p.index;
    }

    return result;
}


/*
 * Params:
 *  initial_angular :   the angular velocity with which the toppling domino
 *      chain was initiated
 *  &lambdas        :   vector with spacings for which to calculate the velocities
 *  mu              :   double, coefficient of friction
 *  full_output :   bool, whether to retreive additional information from
 *      integrator and store it in m_full_output_vec
 *
 * Returns:
 *  double_vec_2d result,
 *   result[0][i]   :   spacial coordinate x
 *   result[0][i]   :   angular velocity of the i-th domino
 *   result[1][i]   :   transversal velocity at the i-th domino
 *
 */
const double_vec_2d DominoChain::make_velocity_array(
        const double initial_angular,
        const double_vec& lambdas,
        const double mu,
        const bool full_output,
        const bool times_only )
{
    if ( full_output )
        m_full_output_vec.clear();

    double_vec_2d result( times_only ? 1 : 3 );
    if ( !times_only )
        result[0].push_back( 0 );   // initial x coordinate

    double psi;
    params p;
    p.index = 0;
    p.angular = initial_angular;
    GslQuad< double_func > integrator( theta_dot_wrapper, m_limit );

    for ( double lambda : lambdas )
    {
        psi = ( p.index == lambdas.size()-1 ) ? M_PI/2 : _psi( lambda );
        p.eta = _eta( lambda );

        if ( times_only )
            result[0].push_back(
                    integrator.integrate( p, 0, psi,
                        m_epsabs, m_epsrel, full_output ));
        else
        {
            if ( p.index != 0 )     // skip first, already added
                result[0].push_back( result[0].back() + lambda + m_h );

            result[1].push_back( p.angular );
            result[2].push_back(
                    transversal( integrator, p, lambda, psi,
                        full_output ));
        }

        p.angular = angular_next( p.index, p.angular, p.eta,
                _theta_hat(lambda), _R(lambda, mu) );
        ++ p.index;
    }

    return result;
}


/*
 * Wrapper of private member function for easy use by user
 */
double DominoChain::intrinsic_angular(
        const double lambda,
        const double mu ) const
{
    return intrinsic_angular(
            _eta( lambda ),
            _theta_hat( lambda ),
            _R( lambda, mu ));
}


double DominoChain::intrinsic_transversal(
        const double lambda,
        const double mu,
        const bool full_output,
        const bool times_only )
{
    params p;
    p.index = m_N;
    p.eta = _eta( lambda );
    p.angular = intrinsic_angular(
            p.eta,
            _theta_hat( lambda ),
            _R( lambda, mu ));
    GslQuad<double_func> integrator( theta_dot_wrapper, m_limit );
    const double time = integrator.integrate(
            p, 0, _psi( lambda ), m_epsabs, m_epsrel, full_output );

    if ( times_only )
        return time;

    if ( full_output )
    {
        m_full_output_vec.clear();
        m_full_output_vec.push_back( integrator.get_result() );
    }

    return ( lambda + m_h ) / time;
}


#if VIDEO

int DominoChain::make_video(
        const std::string filename,
        const double initial_angular,
        const double lambda,
        const double mu,
        const int number_of_pieces,
        const double fps,
        const int length,
        const int width ) const
{
    if ( length % number_of_pieces != 0 )
    {
        std::cerr << "ERROR: `length` needs to be divisible by"
            << "`number_of_pieces` without remainder!\n";
        return -1;
    }

    const double_vec times = _get_times_between_collisions(
            initial_angular, lambda, number_of_pieces, mu );

    cv::VideoWriter writer;
    if ( _open_writer(writer, filename, fps, cv::Size(length, width)) )
        return -1;

    DPRINT( length << " x " << width );

    const int pixelwidth_per_piece = length / number_of_pieces;
    double d_theta;
    double psi = _psi( lambda );
    const double xi_hat_rel = _xi( _theta_hat( lambda ) ) / m_L;
    bool nan_occured = false;
    for ( int index = 0; index < number_of_pieces; ++index )
    {
        if ( gsl_isnan(times[index]) )
        {
            nan_occured = true;
            writer << cv::Mat::ones( width, length, CV_64F /*double*/ );
            continue;
        }

        if ( index == number_of_pieces-1)
            psi = M_PI/2;

        d_theta = psi / ( fps * times[index] );
        for ( double theta = 0; theta < psi; theta += d_theta )
            writer << _make_frame(
                    theta, length, width, index,
                    pixelwidth_per_piece, _eta(lambda), xi_hat_rel );
    }
    if ( nan_occured )
        std::cerr << "WARNING: NaNs occured during calculation!\n";
    std::cerr << "Finished writing video: " << filename << "\n";

    // writer is closed automatically when going out of scope
    return 0;
}


int DominoChain::make_video(
        const std::string filename,
        const double initial_angular,
        const double_vec& lambdas,
        const double mu,
        const double fps,
        const int length,
        const int width ) const
{
    const int number_of_pieces = lambdas.size();
    if ( length % number_of_pieces != 0 )
    {
        std::cerr << "ERROR: `length` needs to be divisible by"
            << "size of `lambdas` without remainder!\n";
        return -1;
    }
    const double_vec times = _get_times_between_collisions(
            initial_angular, lambdas, mu );

    cv::VideoWriter writer;
    if ( _open_writer(writer, filename, fps, cv::Size(length, width)) )
        return -1;

    const int pixelwidth_per_piece = length / number_of_pieces;
    double d_theta;
    double psi;
    bool nan_occured = false;
    for ( int index = 0; index < number_of_pieces; ++index )
    {
        if ( gsl_isnan(times[index]) )
        {
            nan_occured = true;
            writer << cv::Mat::ones( width, length, CV_64F /*double*/ );
            continue;
        }
        psi = ( index == number_of_pieces-1 ) ? M_PI/2 : _psi( lambdas[index] );
        d_theta = psi / ( fps * times[index] );
        for ( double theta = 0; theta < psi; theta += d_theta )
            writer << _make_frame(
                    theta, length, width, index,
                    pixelwidth_per_piece, 0, 0, &lambdas );
    }
    if ( nan_occured )
        std::cerr << "WARNING: NaNs occured during calculation!\n";
    std::cerr << "Finished writing video: " << filename << "\n";

    // writer is closed automatically when going out of scope
    return 0;
}

#endif


std::vector<result>& DominoChain::get_full_output( void )
{
    return m_full_output_vec;
}

void DominoChain::set_pieces_to_be_considered( const int value )
{
    m_N = value;
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
    const double k_value = k( index, theta, eta );
    const double P_over_K_value = P_over_K( theta, angular );
    const double theta_dot = angular * std::sqrt(
            (k_value / (k_value - 1))
            * (1 - P_over_K_value / k_value) );
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

double DominoChain::intrinsic_angular(
        const double eta,
        const double theta_hat,
        const double R) const
{
    const double k_value = k( m_N, 0.0, eta );
    const double return_val = m_omega * std::sqrt( ((k_value-1)/k_value)
            * ( 2 * ( std::cos(m_phi) - std::cos(theta_hat - m_phi)) )
            / (k_value * R * R - k_value + 1) );
	return return_val;
}


/*
 * Transversal velocity V = (λ+h)/∫dθ/dθ'
 */
double DominoChain::transversal(
        GslQuad<double_func>& integrator,
        const params& p,
        const double lambda,
        const double psi,
        const bool full_output )
{
    const double result = ( lambda + m_h ) / integrator.integrate( p, 0, psi,
            m_epsabs, m_epsrel, full_output );

    if ( full_output )
        m_full_output_vec.push_back( integrator.get_result() );

    return result;
}


double DominoChain::angular_next(
        int index,
        double initial_val,
        double eta,
        double theta_hat,
        double R ) const
{
    const double k_value = k( index, 0.0, eta );
    const double P_over_K_value = P_over_K( theta_hat, initial_val );
    return initial_val * std::sqrt( ((k_value-1)/k_value)
            * (1 - P_over_K_value/k_value) ) / R;
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
    const int pieces_to_consider = ( piece_index < m_N ) ? piece_index+1 : m_N;

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
#if WITH_R
    const double xi = _xi( _psi(lambda) );
    return 1 + ( xi + mu * lambda ) / ( xi - mu * m_h );
#else
    DPRINT( "R = 1" );
    return 1;
#endif
}

double DominoChain::_theta_hat( const double lambda ) const
{
    return std::acos( m_h / (m_h + lambda) );
}

double DominoChain::_eta( const double lambda ) const
{
    return ( lambda + m_h ) / m_L;
}


#if VIDEO

// video helper methods
// static
int DominoChain::_open_writer(
        cv::VideoWriter& writer,
        const std::string filename,
        const double fps,
        const cv::Size framesize )
{
    writer.open(
            filename,
            // CV_FOURCC('D', 'I', 'V', 'X'),
            -1,
            fps,
            framesize,
            false /*isColor*/ );
    if ( !writer.isOpened() )
    {
        std::cerr << "ERROR: Failed to open `VideoWriter`\n";
        return -1;
    }
    return 0;
}

// non-static
const cv::Mat DominoChain::_make_frame(
        double theta,
        const int length,
        const int width,
        const int index,
        const int pixelwidth_per_piece,
        const double eta,
        const double min_height,
        const double_vec* lambdas ) const
{
    double_vec heights( length * width, 1 );
    const int step = width * pixelwidth_per_piece;
    double xi_rel;
    double eta_new;

    for ( int i = index; i >= 0; --i )
    {
        xi_rel = _xi( theta ) / m_L;
        // if `min_height` is not 0, check if there is no more significant
        // change in domino positions. In that case assign the remaining pieces
        // with the corresponding value and exit the iteration
        if ( min_height && std::fabs( xi_rel - min_height ) < 0.004 /*≈1/256*/ )
        {
            std::fill_n( heights.begin(), (i+1) * step, min_height );
            break;
        }
        std::fill_n( heights.begin() += (i * step ), step, xi_rel );
        // θ_i+1 = arcsin(η * cos(θ_i) - h/L) + θ_i
        eta_new = lambdas ? _eta( lambdas->at(i) ) : eta;
        theta = std::asin( eta_new * std::cos( theta ) - m_h / m_L ) + theta;
    }

    cv::Mat frame( heights );
    frame = frame.reshape( 0, length );

    return frame.t();
}

const double_vec DominoChain::_get_times_between_collisions(
        const double initial_angular,
        const double lambda,
        const int length,
        const double mu ) const
{
    return const_cast<DominoChain*>(this)->make_velocity_array(
            initial_angular,
            lambda,
            length,
            mu,
            false /*full_output*/,
            true /*times_only*/ )[0];
}

const double_vec DominoChain::_get_times_between_collisions(
        const double initial_angular,
        const double_vec& lambdas,
        const double mu ) const
{
    return const_cast<DominoChain*>(this)->make_velocity_array(
            initial_angular,
            lambdas,
            mu,
            false /*full_output*/,
            true /*times_only*/ )[0];
}

#endif

