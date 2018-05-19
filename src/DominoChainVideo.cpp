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
    cv::VideoWriter writer(
            filename,
            CV_FOURCC('D', 'I', 'V', 'X'),
            m_fps,
            m_size,
            false   // isColor
            );
    // TODO
    // What about `writer.open()`?
    // check `writer.isOpen()`
    double d_theta;
    double psi = _psi( lambda );
    for ( int index = 0; index < m_length; ++index )
    {
        d_theta = psi / ( m_fps * times[index] );
        for ( double theta = 0; theta < psi; theta += d_theta )
            writer << _make_frame( index, theta );
    }

    // writer is closed automatically when going out of scope
    // TODO
    // What about errors?
    return 0;
}


//////////////////////
// Private Members  //
//////////////////////

cv::Mat DominoChainVideo::_make_frame(
        const int index,
        const double theta ) const
{
    // double heights[m_length];
    auto heights = std::make_unique<double[]>(m_length);
    memset( &heights, 1, m_length );
    for ( int j = index; j >= 0; ++j )
    {
        // TODO
        // check if |ξ_i - ξ^| < 0.004 (≈1/256)
        heights[j] = _xi( theta ) / m_L;
    }

    auto mat_data = std::make_unique<double[]>(m_length*m_width);
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
