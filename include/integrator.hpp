#pragma once

#include "gsl_function_cpp.hpp"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <memory>
#include <functional>

typedef std::unique_ptr< gsl_integration_workspace,
        std::function<void(gsl_integration_workspace*)>
        > gsl_integration_workspace_cpp;
typedef std::unique_ptr<gsl_function> gsl_func_ptr;

/*
 * wrapper for gsl_integration
 *  https://scicomp.stackexchange.com/a/27248
 *  https://stackoverflow.com/a/24151084/9133910
 *
 */
template<typename F>
class gsl_quad
{
    public:
        gsl_quad(F func, size_t limit)
            : m_func(func), m_limit(limit),
            m_workspace(gsl_integration_workspace_alloc(limit),
                    gsl_integration_workspace_free)
        {}

        double integrate(double min, double max, double epsabs, double epsrel)
        {
            auto func_ptr = [this](double x){return this->m_func(x);};
            gsl_function_cpp<decltype(func_ptr)> f(func_ptr);
            // unique_ptr?
            gsl_func_ptr gsl_f = static_cast<gsl_func_ptr>(&f);
        }

    private:
        F m_func;
        size_t m_limit;
        gsl_integration_workspace_cpp m_workspace;
};
