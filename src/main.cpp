#include "../include/DominoChain.hpp"
#include <iostream>
#include <assert.h>


int main()
{
    domino d;
    d.height = 0.02;    // m
    d.width = 0.004;    // m

    double_vec_2d intrinsic_velocities;
    DominoChain dc (d, 10, 100, 1000);

    double_vec lambdas {0.004, 0.005, 0.006, 0.007, 0.008, 0.009,
		0.01, 0.011};

    try
    {
        intrinsic_velocities = dc.make_velocity_array(lambdas, 0.2);
    }
    catch ( std::runtime_error &err)
    {
        std::cerr << err.what() << std::endl;
    }
    assert (intrinsic_velocities[0].size() == intrinsic_velocities[1].size());

    std::cout << "phi_dot\tV\n";
    for ( size_t i = 0; i < intrinsic_velocities[0].size(); ++i)
        std::cout << intrinsic_velocities[0][i] << "\t"
            << intrinsic_velocities[1][i] << "\n";

    std::cout << "\n\n" << std::endl;

    double_vec_2d velocity_at_x;
    try
    {
        velocity_at_x = dc.make_velocity_array(2, 0.008, 0.2);
    }
    catch ( std::runtime_error &err )
    {
        std::cerr << err.what() << std::endl;
    }
    assert ( velocity_at_x[0] == velocity_at_x[1] );
    assert ( velocity_at_x[1] == velocity_at_x[2] );

    std::cout << "phi_dot\tV\n";
    for ( size_t i = 0; i < velocity_at_x[0].size(); ++i )
        std::cout << velocity_at_x[0][i] << "\t"
            << velocity_at_x[1][i] << "\t"
            << velocity_at_x[2][i] << "\n";

    std::cout << std::endl;

    return 0;
}
