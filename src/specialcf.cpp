#include <fem.hpp>

namespace f2c {
#include <f2c.h>
}

extern "C" /* Subroutine */ int zbesi_(f2c::doublereal *zr, f2c::doublereal *zi, const f2c::doublereal *fnu, const f2c::integer *kode, const f2c::integer *n, f2c::doublereal *cyr, f2c::doublereal *cyi, f2c::integer * nz, f2c::integer *ierr);
extern "C" /* Subroutine */ f2c::doublereal gamln_(f2c::real *z__, f2c::integer *ierr);

Complex zbesi( Complex z, double order, int kode_ )
{
  // input arguments 
  f2c::doublereal fnu = order;
  f2c::integer kode = kode_;
  f2c::integer n = 1;
  f2c::doublereal zr = z.real();
  f2c::doublereal zi = z.imag();

  // outputs
  f2c::doublereal cyr, cyi;
  f2c::integer nz;
  f2c::integer ierr;

  zbesi_(&zr, &zi, &fnu, &kode, &n, &cyr, &cyi, &nz, &ierr);
  if(nz>0) cout << "Number of underflows: " << nz << endl;
  if(ierr>0) cout << "Error: " << ierr << endl;

  return Complex(cyr, cyi);
}

double gamln(double x_)
{
  f2c::real x = x_;
  f2c::integer ierr;
  double res =  gamln_(&x, &ierr);
  if(ierr>0) cout << "Error: " << ierr << endl;
  return res;
}

namespace ngfem {

    class SpecialCoefficientFunction_ZBESI : public CoefficientFunction
    {
        shared_ptr<CoefficientFunction> arg;

        double order;
        int kode;

      public:

        SpecialCoefficientFunction_ZBESI(shared_ptr<CoefficientFunction> arg_, double order_, int kode_)
          : CoefficientFunction(1, true),
          arg(arg_), order(order_), kode(kode_)
        {
        }

        double Evaluate (const BaseMappedIntegrationPoint & ip) const override 
        {
            throw Exception("SpecialCoefficientFunction_ZBESI::Evaluate called");
        }


        virtual void Evaluate(const BaseMappedIntegrationPoint & ip, FlatVector<Complex> result) const override
        {
            Complex z = arg->EvaluateComplex(ip);
            result(0) = zbesi(z, order, kode);
        }
    };




    class SpecialCoefficientFunction_DGAMLN : public CoefficientFunction
    {
        shared_ptr<CoefficientFunction> arg;

      public:
        SpecialCoefficientFunction_DGAMLN(shared_ptr<CoefficientFunction> arg_)
          : CoefficientFunction(1, false),
          arg(arg_)
        { }

        virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override {
            double x = arg->Evaluate(ip);
            return gamln(x);
        }

    };
}

#undef min
#undef max
#undef abs
#include <python_ngstd.hpp>
PYBIND11_MODULE(special_functions, m) {
    m.def("Bessel", [] (shared_ptr<ngfem::CoefficientFunction> arg, double order, int kode) -> shared_ptr<ngfem::CoefficientFunction>
          {
          return make_shared<ngfem::SpecialCoefficientFunction_ZBESI>(arg, order, kode);
          }, py::arg("z"), py::arg("order")=0, py::arg("kode")=1, py::doc(R"DOCSTRING_(

Input
  z      - Complex argument
  order  - DOUBLE PRECISION order>=0
  kode   - A parameter to indicate the scaling option
           kode=1  returns
                   CY=I(order,z)
               =2  returns
                   CY=exp(-abs(X))*I(order,z)
                   where X=Re(z)
)DOCSTRING_")
    );
    m.def("Gamma", [] (shared_ptr<ngfem::CoefficientFunction> arg) -> shared_ptr<ngfem::CoefficientFunction>
          {
          return make_shared<ngfem::SpecialCoefficientFunction_DGAMLN>(arg);
          });
}
