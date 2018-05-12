#include "../include/GslQuad.hpp"
#include <iostream>
#include <cmath>

typedef std::pair<double, double> tuple;
typedef struct params_struct
{
    double c;
} params;


template<typename F, typename P>
double doit(F func, tuple const& range, P params,
        double epsabs = 1.49e-8, double epsrel = 1.49e-8,
        int limit = 100)
{
    return GslQuad<F>(func, limit).integrate(params,
            range.first, range.second, epsabs, epsrel);
}


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
    auto h = [](double x, params p){return 1/std::pow(x, p.c);};
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
        << "\\int_0^1 dx/x = "
        << doit(h, {0, 1}, p) << '\n'
        << "\\int_{-1}^1 sin^2(x)/x^2 dx = "
        << doit([](double x, params p){return std::sin(x)*std::sin(x)/(x*x);}, {-1, 1}, p) << '\n'
        << std::endl;
    }
    catch (std::runtime_error& err)
    {
        std::cout << err.what() << std::endl;
    }

    return 0;
}
