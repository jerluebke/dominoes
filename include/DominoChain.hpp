#pragma once

#ifndef VIDEO
#define VIDEO 0
#endif

#include "GslQuad.hpp"
#include <vector>
#include <string>

#if VIDEO
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/videoio.hpp>
#endif


typedef std::vector<double> double_vec;
typedef std::vector< std::vector<double> > double_vec_2d;

typedef struct domino_struct
{
    double height;
    double width;
} domino;

typedef struct params_struct
{
    int index;
    double eta;
    double angular;
} params;

typedef std::function< double (double, params) > double_func;


class DominoChain
{
    public:
        DominoChain (
                const domino& d,
                int N = 10,
                const int limit = 100,
                const double epsabs = 1.49e-8,
                const double epsrel = 1.49e-8 );

        // velocity vectors
        const double_vec_2d make_velocity_array(
                const double_vec& lambdas,
                const double mu,
                const bool full_output = false,
                const bool times_only = false );

        const double_vec_2d make_velocity_array(
                const double initial_angular,
                const double lambda,
                const int number_of_pieces,
                const double mu,
                const bool full_output = false,
                const bool times_only = false );

        const double_vec_2d make_velocity_array(
                const double initial_angular,
                const double_vec& lambdas,
                const double mu,
                const bool full_output = false,
                const bool times_only = false );

        // for fitting μ
        double intrinsic_angular(
                const double lambda,
                const double mu ) const;

        double intrinsic_transversal(
                const double lambda,
                const double mu,
                const bool full_output = false,
                const bool times_only = false );

#if VIDEO
        int make_video(
                const std::string filename,
                const double initial_angular,
                const double lambda,
                const double mu,
                const int number_of_pieces = 128,
                const double fps = 30,
                const int length = 512,
                const int width = 64 ) const;

        int make_video(
                const std::string filename,
                const double initial_angular,
                const double_vec& lambdas,
                const double mu,
                const double fps = 30,
                const int length = 512,
                const int width = 64 ) const;
#endif

        std::vector<result>& get_full_output( void );

		void set_pieces_to_be_considered( const int value );


    private:
        // function to be integrated
        double theta_dot(
                const double theta,
                const int index,
                const double eta,
                const double angular ) const;

        const double_func theta_dot_wrapper;

        // full output of integration
        std::vector<result> m_full_output_vec;

        // integrator params
        const size_t m_limit;
        const double m_epsabs;
        const double m_epsrel;

        // domino quantities
        int m_N;          // number of dominoes to be considered
        const double m_L;       // domino height
        const double m_h;       // domino width
        const double m_phi;     // = arctan(h/L)
        const double m_omega;   // eigenfrequency

        // used to initialize m_phi and m_omega
        static double _phi( const double L, const double h );
        static double _omega( const double L, const double phi );

        // integration helper methods
        double intrinsic_angular(
                const double eta,
                const double theta_hat,
                const double R ) const;

        double transversal(
                GslQuad<double_func>& integrator,
                const params& p,
                const double lambda,
                const double psi,
                const bool full_output = false );

        double angular_next(
                const int index,
                const double initial_val,
                const double eta,
                const double theta_hat,
                const double R ) const;

        double P_over_K(
                const double theta,
                const double intitial_angular_speed ) const;

        double k(
                const int piece_index,
                const double theta_initial,
                const double eta ) const;

        double theta_dot_rel(
                const double theta,
                const double eta ) const;

        // instance helper methods
        double _psi( const double lambda ) const;   // angle of impact
        double _xi( const double psi ) const;       // height of impact
        double _R( const double lambda,
                const double mu ) const;            // transmission of angular speed
        double _theta_hat( const double lambda ) const;   // final angle
        double _eta( const double lambda ) const;   // = (λ+h)/L

#if VIDEO
        // video helper methods
        static int _open_writer(
                cv::VideoWriter& writer,
                const std::string filename,
                const double fps,
                const cv::Size framesize );

        const cv::Mat _make_frame(
                double theta,
                const int length,
                const int width,
                const int index,
                const int pixelwidth_per_piece,
                const double eta,
                const double min_height = 0,
                const double_vec* lambdas = nullptr ) const;

        const double_vec _get_times_between_collisions(
                const double initial_angular,
                const double lambda,
                const int length,
                const double mu ) const;

        const double_vec _get_times_between_collisions(
                const double initial_angular,
                const double_vec& lambdas,
                const double mu ) const;
#endif

};
