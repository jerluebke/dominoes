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
class GslQuad
{
    public:
        GslQuad(F func, size_t limit)
            : m_gsl_wrapper(func)
            , m_limit(limit)
            , m_workspace(gsl_integration_workspace_alloc(limit),
                    gsl_integration_workspace_free)
        {
            // turn error handler off to prevent termination of program
            gsl_set_error_handler_off();
        }

        double integrate(double min, double max, double epsabs, double epsrel)
        {
            // cast c++ wrapper to gsl_function
            gsl_function* gsl_f = static_cast<gsl_function*>(&m_gsl_wrapper);

            // do integration
            double result, error;
            int status = gsl_integration_qags( gsl_f, min, max,
                    epsabs, epsrel, m_limit, m_workspace.get(),
                    &result, &error );

            if (status)
                _handle_error(status);

            return result;
        }

    private:
        size_t m_limit;
        GslFunctionCpp<F> m_gsl_wrapper;
        gsl_integration_workspace_cpp m_workspace;

        void _handle_error(int status) const
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
    return GslQuad<F>(func, limit).integrate(
            range.first, range.second, epsabs, epsrel);
}
