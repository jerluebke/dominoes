#include "../include/DominoChainVideo.hpp"
#include <array>
#include <memory>


//////////////////
// Constructor  //
//////////////////

DominoChainVideo::DominoChainVideo(
        const domino& d,
        double fps,
        const int length,
        const int width,
        int N )
    : DominoChain( d, N )
    , m_fps( fps )
    , m_length( length )
    , m_width( width )
    , m_size( width, length )
    {}


//////////////////////
// Public Members   //
//////////////////////

int DominoChainVideo::make_video(
        const std::string filename,
        const double initial_angular,
        const double lambda,
        const double mu ) const
{
    double_vec times = _get_times_between_collisions(
            initial_angular, lambda, mu );
    cv::VideoWriter writer;
    writer.open(
            filename,
            CV_FOURCC('D', 'I', 'V', 'X'),
            m_fps,
            m_size,
            false /*isColor*/ );
    if ( !writer.isOpened() )
    {
        std::cerr << "ERROR: Failed to open `VideoWriter`\n";
        return -1;
    }
    double d_theta;
    double psi = _psi( lambda );
    double xi_hat_rel = _xi( _theta_hat( lambda ) ) / m_L;
    for ( int index = 0; index < m_length; ++index )
    {
        d_theta = psi / ( m_fps * times[index] );
        for ( double theta = 0; theta < psi; theta += d_theta )
            writer << _make_frame( index, theta, xi_hat_rel );
    }

    // writer is closed automatically when going out of scope
    return 0;
}


//////////////////////
// Private Members  //
//////////////////////

cv::Mat DominoChainVideo::_make_frame(
        const int index,
        const double theta,
        const double min_height ) const
{
    auto heights = std::make_unique<double[]>( m_length );
    memset( &heights, 1, m_length );

    double xi_rel;
    for ( int j = index; j >= 0; --j )
    {
        xi_rel = _xi( theta ) / m_L;
        if ( std::fabs( xi_rel - min_height ) < 0.004 /*≈1/256*/ )
        {
            memset( &heights[0], min_height, index );
            break;
        }
        else
            heights[j] = _xi( theta ) / m_L;
    }

    auto mat_data = std::make_unique<double[]>( m_length*m_width );
    memset( &mat_data, 1, m_length*m_width );

    double start = 0;
    for ( int k = 0; k <= index; ++k )
    {
        memset( &mat_data[start], heights[k], m_width );
        start += m_width;
    }

    cv::Mat frame ((std::vector<double>( *mat_data.get() )));
    frame.reshape( 0, m_width );

    return frame.t();
}


double_vec DominoChainVideo::_get_times_between_collisions(
        const double initial_angular,
        const double lambda,
        const double mu ) const
{
    double_vec_2d result = make_velocity_array(
            initial_angular,
            lambda,
            m_length,
            mu,
            false /*full_output*/,
            true /*times_only*/ );
    return result[0];
}
