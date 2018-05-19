#include "../include/DominoChain.hpp"

#include <opencv2/core.hpp>
#include <opencv2/videoio.hpp>
#include <opencv2/highgui.hpp>

#include <iostream>
#include <vector>
#include <cmath>


std::vector<double> make_content(
        int size,
        double psi,
        double d_theta,
        int index,
        double eta,
        double h,
        double L,
        DominoChain dc )
{
    std::vector<double> result ( size, 1 );
    double theta_i;

    for ( double theta_0 = 0; theta_0 < psi; theta_0 += d_theta )
    {
        theta_i = theta_0;
        for ( int j = index; j >= 0; --j )
        {
            result[j] = dc._xi( theta_i );
            theta_i = std::asin( eta * std::cos( theta_i ) - h/L ) + theta_i;
        }
    }
}


cv::Mat make_frame( double d_theta, int index, cv::Size framesize )
{
    // std::vector<double> content = make_content(  );
}


int main()
{
    domino d;
    d.height = 2;
    d.width = 0.4;
    DominoChain dc ( d );
    double initial_velocity = 6.0;
    double spacing = 1.0;
    double mu = 0.2;
    int number_of_pieces = 512;
    // TODO: only set first 10 elements, fill rest with intrinsic speed
    // implement in that in `DominoChain` as optimization
    double_vec times = dc.get_times_between_collisions(
            initial_velocity,
            spacing,
            number_of_pieces,
            mu );
    double psi = dc._psi( spacing );

    double fps = 60;
    int codec = CV_FOURCC('D', 'I', 'V', 'X');
    cv::String filename = "./test_video.avi";
    cv::Size framesize( number_of_pieces, 100 );
    cv::VideoWriter writer;
    writer.open( filename, codec, fps, framesize, false );

    if ( !writer.isOpened() )
    {
        std::cerr << "Failed to open video writer!\n";
        return -1;
    }

    for ( int i = 0; i < number_of_pieces; ++i )
        writer << make_frame( psi / ( fps * times[i] ), i, framesize );

    return 0;
}
