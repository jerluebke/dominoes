#include "../include/DominoChain.hpp"

/*
 * TODO:
 *  testing and optimizing, perhaps multithreading
 *
 *  make video
 *
 *  add least-square-fitting (which algorithm is used by numpy/scipy?)
 *
 */

int main()
{
    domino d;
    d.height = 0.02;    // m
    d.width = 0.004;    // m
    double_vec mus { 0.1, 0.2, 0.4, 0.6 };
    DominoChain dc (d, 10, 20, 100);
    double_vec_2d intrinsic_velocities;

    double_vec lambdas {
        // 0.0,
        0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009,
        0.01, 0.011, 0.012, 0.013, 0.014, 0.015, 0.016, 0.017, 0.018, 0.019
        // , 0.02
    };

    std::cout << "intrinsic velocities:\n";

    for ( double mu : mus )
    {
        std::cout << "\\mu = " << mu << "\n";

        intrinsic_velocities = dc.make_velocity_array(lambdas, mu, false);
        std::cout << "l\tphi_dot\tV\n";
        for ( size_t i = 0; i < lambdas.size(); ++i)
            std::cout << lambdas[i] << '\t'
                << intrinsic_velocities[0][i] << "\t"
                << intrinsic_velocities[1][i] << "\n";
                // << dc.get_full_output( i ).str();

        std::cout << "\n\n";
    }

    double_vec initial_velocities { 2*M_PI/3, M_PI, 3*M_PI/2, 2*M_PI, 5*M_PI/2 };
    double_vec_2d velocity_at_x;

    std::cout << "velocity at x:\n"
        << "\\mu = " << 0.2 << ", \\lambda = " << 0.006 << "\n";

    for ( double initial : initial_velocities )
    {
        std::cout << "initial angular velocity = " << initial << "\n";

        velocity_at_x = dc.make_velocity_array(initial, 0.006, 20, 0.2, true);
        std::cout << "x\tphi_dot\tV\n";
        for ( size_t i = 0; i < velocity_at_x[0].size(); ++i )
            std::cout << velocity_at_x[0][i] << "\t"
                << velocity_at_x[1][i] << "\t"
                << velocity_at_x[2][i] << "\n";
                // << dc.get_full_output( i ).str();

        std::cout << "\n\n";
    }

    return 0;
}
