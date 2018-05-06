# Integrator
---
### Core
GSL QAGS (adaptive integration with singularities)  

    int gsl_integration_qags(const gsl_function* f,
                             double a,
                             double b,
                             double epsabs,
                             double epsrel,
                             size_t limit,
                             gsl_integration_workspace* w,
                             double* result,
                             double* abserr)

workspace -> use `unique_ptr`

    gsl_integration_workspace giw
    giw* gsl_integration_workspace_alloc(size_t n)
    void gsl_integration_workspace_free(giw* w)

function

    gsl_function gf
    double (*function) gf.function(double x, void* params)
    void* gf.params
