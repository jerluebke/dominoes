#include "./include/FunctionTimer.h"
#include <vector>
#include <cmath>
#include <iostream>

typedef std::vector<std::vector<double>> double_vec_2d;

double_vec_2d make_3xN_vec( const int N )
{
    double_vec_2d result;
    result.resize( 3 );
    for ( int i = 0; i < N; ++i)
    {
        result[0].push_back( std::sin(i) );
        result[1].push_back( std::cos(i) );
        result[2].push_back( std::tan(i) );
    }
    return result;
}

double_vec_2d make_Nx3_vec( const int N )
{
    double_vec_2d result;
    result.resize( N );
    for ( int i = 0; i < N; ++i )
    {
        result.push_back({
                std::sin(i),
                std::cos(i),
                std::tan(i),
                });
    }
    return result;
}

int main()
{
    int N = 100;
    double time_first;
    {
        BlockTimer timer( time_first );
        double_vec_2d res_first = make_3xN_vec( N );
    }
    std::cout << "3x" << N << " : " << time_first << "ms\n";

    double time_second;
    {
        BlockTimer timer( time_second );
        double_vec_2d res_second = make_Nx3_vec( N );
    }
    std::cout << N << "x3 : " << time_second << "ms\n";

    return 0;
}
