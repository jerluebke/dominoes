#pragma once

#include <gsl/gsl_math.h>

/*
 * wrapper for gsl_function to use with member function of some class
 *  https://stackoverflow.com/a/18181494/9133910
 *
 * usage:
 *  SomeClass* ptr1 = this;
 *  auto ptr2 = [=](double x)->double{return ptr1->some_member(x);};
 *  gsl_function_cpp<decltype(ptr2)> f(ptr2);
 *  gsl_function *F = static_cast<gsl_function*>(&f);
 *  // do something with F using gsl routines ...
 *
 * Notes:
 *  use template for performance increase compared to std::function
 *
 */
template<typename F>
class GslFunctionCpp : public gsl_function
{
    public:
        GslFunctionCpp(const F& func)
            : m_func(func)
        {
            function = &GslFunctionCpp::invoke;
            params = this;
        }


    private:
        const F& m_func;
        static double invoke(double x, void* params)
        {
            return static_cast<GslFunctionCpp*>(params)->m_func(x);
        }
};
