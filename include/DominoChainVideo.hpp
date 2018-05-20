#pragma once

#include "DominoChain.hpp"
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/videoio.hpp>
#include <string>


class DominoChainVideo : public DominoChain
{
    public:
        DominoChainVideo(
                const domino& d,
                double fps = 30,
                const int length = 512,
                const int width = 64,
                int N = 10 );

        int make_video(
                const std::string filename,
                const double initial_angular,
                const double lambda,
                const double mu );

        double_vec_2d make_velocity_array(
                const double initial_angular,
                const double lambda,
                const int number_of_pieces,
                const double mu,
                const bool full_output = false,
                const bool times_only = false );

        double m_fps;

    private:
        const int m_length;
        const int m_width;
        const cv::Size m_size;

        cv::Mat _make_frame(
                const int index,
                const double theta,
                const double min_height ) const;

        double_vec _get_times_between_collisions(
                const double initial_angular,
                const double lambda,
                const double mu );

};
