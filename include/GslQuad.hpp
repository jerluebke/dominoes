#pragma once

#include "gsl_function_cpp.hpp"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <memory>
#include <functional>
#include <sstream>
#include <iostream>

typedef std::unique_ptr< gsl_integration_workspace,
        std::function<void(gsl_integration_workspace*)>
        > gsl_integration_workspace_cpp;

typedef struct result
{
    // double integration_result;
    double error;
    int status;
    std::string errormsg;

    std::string str()
    {
        std::stringstream output;
        output << "{ "
            // << "\"result\" : " << integration_result << ",\n"
            << "\"error\" : " << error << ",\n"
            << "\"gsl_errno\" : " << status << ",\n"
            << "\"error_msg\" : \"" << errormsg << "\" }\n"
            << '\n';
        return output.str();
    }

} result;

/*
 * TODO:
 * This doesn't work - Linke Error
 *
 * clang output:
 * DominoChain-cec92d.o : error LNK2005: "class std::basic_ostream<char,struct std::char_traits<char> > & __cdecl operator<<(class std::basic_ostream<char,struct std::char_traits<char> > &,struct result const &)" (??6@YAAEAV?$basic_ostream@DU?$char_traits@D@std@@@std@@AEAV01@AEBUresult@@@Z) already defined in main-bc4c96.o
 *
 */
// std::ostream& operator<< ( std::ostream& os, const result& rs )
// {
//     os << "{ result = " << rs.integration_result << ",\n"
//         << "error = " << rs.error << ",\n"
//         << "gsl errno: " << rs.status << ",\n"
//         << "error msg: " << rs.errormsg << " }\n"
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
                const double epsrel = 1.49e-8,
                const bool full_output = false)
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

            if ( full_output )
            {
                // result_struct.integration_result = result;
                result_struct.error = error;
                result_struct.status = status;
                result_struct.errormsg = std::string( gsl_strerror(status) );
            }

            return result;
        }


        result get_result() const
        { return result_struct; }


    private:
        const F m_func;
        const size_t m_limit;
        gsl_integration_workspace_cpp m_workspace;
        result result_struct;

};

