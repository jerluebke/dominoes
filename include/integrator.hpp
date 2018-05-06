#pragma once

#include "gsl_function_cpp.hpp"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <memory>
#include <functional>
#include <stdexcept>
#include <sstream>
#include <string>
#include <utility>

typedef std::unique_ptr< gsl_integration_workspace,
        std::function<void(gsl_integration_workspace*)>
        > gsl_integration_workspace_cpp;
typedef std::unique_ptr<gsl_function> gsl_func_ptr;
typedef std::pair<double, double> tuple;

/*
 * wrapper for gsl_integration
 *  https://scicomp.stackexchange.com/a/27248
 *
 * using QAGS adaptive integration with singularites
 *  (with fixed borders)
 *
 * error handling:
 *  https://stackoverflow.com/a/24151084/9133910
 *
 */
template<typename F>
class gsl_quad
{
    public:
        gsl_quad(F func, size_t limit)
            : m_func(func), m_limit(limit)
            , m_workspace(gsl_integration_workspace_alloc(limit),
                    gsl_integration_workspace_free)
        {}

        double integrate(double min, double max, double epsabs, double epsrel)
        {
            // turn error handler off to prevent termination of program
            gsl_set_error_handler_off();

            // bind member m_func to make cpp wrapper
            auto func_ptr = [this](double x)->double
                {return this->m_func(x);};
            gsl_function_cpp<decltype(func_ptr)> f(func_ptr);
            // cast to gsl_function
            gsl_function* gsl_f = static_cast<gsl_function*>(&f);

            // do integration
            // gsl_f is already a pointer
            double result, error;
            int status = gsl_integration_qags( gsl_f, min, max,
                    epsabs, epsrel, m_limit, m_workspace.get(),
                    &result, &error );

            if (status)
                handleError(status);

            return result;
        }

    private:
        F m_func;
        size_t m_limit;
        gsl_integration_workspace_cpp m_workspace;

        void handleError(int status) const
        {
            std::stringstream msg;
            msg << "GSL ERROR: " << std::string(gsl_strerror(status));
			msg << " (errno " << status << ")";
            throw std::runtime_error(msg.str());
        }

};


template<typename F>
double doit(F func, tuple const& range,
        double epsabs = 1.49e-8, double epsrel = 1.49e-8,
        int limit = 1000)
{
    return gsl_quad<F>(func, limit).integrate(
            range.first, range.second, epsabs, epsrel);
}
