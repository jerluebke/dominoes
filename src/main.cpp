#include "../include/integrator.hpp"
#include <iostream>
#include <cmath>

typedef struct params_struct
{
    double c;
} params;

int main()
{
    // auto f = [](double x)->double{return -(x*x)+1;};
    // GslQuad<decltype(f)> integrator(f, 100);
    params p;
    p.c = 2;
    params q;
    q.c = 5;
    auto f = [](double x, params p){return -(x*x)+p.c;};
    auto g = [](double x, params p){return std::exp(-(x*x)/p.c);};
    GslQuad<decltype(g)> integrator(g, 100);
    try
    {
    std::cout
        << "\\int_-1^1 -x^2+1 dx = "
        << doit<decltype(f), decltype(p)>(f, {-1, 1}, p) << '\n'
        // << integrator.integrate(-1, 1) << '\n';
        << "\\int_-1^1 \\exp(-x^2/2) dx = "
        << integrator.integrate<params>(p, -1, 1) << '\n'
        << "\\int_-1^1 \\exp(-x^2/5) dx = "
        << integrator.integrate<params>(q, -1, 1) << '\n'
        << std::endl;
    }
    catch (std::runtime_error& err)
    {
        std::cout << err.what() << std::endl;
    }
}
