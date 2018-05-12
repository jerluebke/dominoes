#pragma once

#include "gsl_function_cpp.hpp"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <memory>
#include <functional>
#include <stdexcept>
#include <sstream>

typedef std::unique_ptr< gsl_integration_workspace,
        std::function<void(gsl_integration_workspace*)>
        > gsl_integration_workspace_cpp;

// typedef struct result
// {
//     double result;
//     double error;
//     int status;
//     std::string errormsg;
// 
// } result;
// 
// std::ostream& operator<< ( std::ostream& os, const result& rs )
// {
//     os << "result = " << rs.result
//         << "error = " << rs.error
//         << "gsl errno: " << rs.status
//         << "error msg: " << rs.errormsg
//         << '\n';
//     return os;
// }


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
        GslQuad( const F func, const size_t limit = 100 )
            : m_func( func )
            , m_limit( limit )
            , m_workspace( gsl_integration_workspace_alloc(limit),
                    gsl_integration_workspace_free )
        {
            // turn error handler off to prevent termination of program
            gsl_set_error_handler_off();
        }

        template<typename P>
        double integrate(
                const P& params,
                const double min, const double max,
                const double epsabs = 1.49e-8,
                const double epsrel = 1.49e-8
                // const bool full_output = false )
               ) const
        {
            // bind member method for c++ wrapper
            auto func_ptr = [this, &params](double x)->double
                { return this->m_func(x, params); };
            // make c++ wrapper
            GslFunctionCpp<decltype(func_ptr)> gsl_wrapper(func_ptr);
            // cast c++ wrapper to gsl_function
            gsl_function* gsl_f = static_cast<gsl_function*>(&gsl_wrapper);

            // do integration
            double result, error;
            int status = gsl_integration_qags( gsl_f, min, max,
                    epsabs, epsrel, m_limit, m_workspace.get(),
                    &result, &error );

            if (status)
                handle_error(status, result, error);
            
            // if ( full_output )
            // {
            //     result_struct.result = result;
            //     result_struct.error = error;
            //     result_struct.status = status;
            //     result_struct.errormsg = std::string( gsl_strerror(status) );
            // }

            return result;
        }


        // result get_result() const
        // { return result_struct; }


    private:
        const F m_func;
        const size_t m_limit;
        gsl_integration_workspace_cpp m_workspace;
        // result result_struct;

        void handle_error( int status, double result, double error ) const
        {
            std::stringstream msg;
            msg << "GSL ERROR: " << std::string(gsl_strerror(status));
			msg << " (errno " << status << ")";
			msg << " - result = " << result << ", error = " << error;
            throw std::runtime_error(msg.str());
        }

};

