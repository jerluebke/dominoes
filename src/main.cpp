#include "../include/DominoChain.hpp"

/* 
 * TODO:
 *  find out valid domain for parameters (hopefully from the mathematics ...)
 *  handle input accordingly (raise exception)
 *
 *  compile and wrap in cython
 *
 *  testing and optimizing, perhaps multithreading
 *
 *  add least-square-fitting (which algorithm is used by numpy/scipy?)
 *
 * */

int main()
{
    domino d;
    d.height = 0.02;    // m
    d.width = 0.004;    // m

    double_vec_2d intrinsic_velocities;
    DominoChain dc (d, 10, 20, 100);

    double_vec lambdas{
        0.0,
        0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009,
        0.01, 0.011, 0.012, 0.013, 0.014, 0.015, 0.016, 0.017, 0.018, 0.019
        , 0.02 
    };

    intrinsic_velocities = dc.make_velocity_array(lambdas, 0.2, true);

    std::cout << "l\tphi_dot\tV\n";
    for ( size_t i = 0; i < lambdas.size(); ++i)
        std::cout << lambdas[i] << '\t'
			<< intrinsic_velocities[0][i] << "\t"
            << intrinsic_velocities[1][i] << "\n"
            << dc.get_full_output( i ).print();

    std::cout << "\n\n";

    double_vec_2d velocity_at_x;

    velocity_at_x = dc.make_velocity_array(2, 0.01, 20, 0.2, true);

    std::cout << "x\tphi_dot\tV\n";
    for ( size_t i = 0; i < velocity_at_x[0].size(); ++i )
        std::cout << velocity_at_x[0][i] << "\t"
            << velocity_at_x[1][i] << "\t"
            << velocity_at_x[2][i] << "\n"
            << dc.get_full_output( i ).print();

    std::cout << std::endl;

    return 0;
}
