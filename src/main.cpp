#include "../include/integrator.hpp"
#include <iostream>
#include <cmath>

int main()
{
    try
    {
    std::cout 
        << "\\int_-1^1 -x^2+1 dx = "
        << doit([](double x){return -(x*x)+1;}, {-1, 1}) << '\n'
        << "\\int_-1^1 \\exp(-x^2/2) dx = "
        << doit([](double x){return std::exp(-(x*x)/2);}, {-1, 1}) << '\n'
        << "\\int_{-\\pi}^{\\pi} \\sin^2(x)/x^2 dx = "
        << doit([](double x){return std::sin(x)*std::sin(x)/(x*x);}, {-M_PI, M_PI}) << '\n'
        << "\\int_0^1 x^{-2} dx = "
        << doit([](double x){return 1/(x*x);}, {0, 1})
        << std::endl;
    }
    catch (std::runtime_error& err)
    {
        std::cout << err.what() << std::endl;
    }
} 
