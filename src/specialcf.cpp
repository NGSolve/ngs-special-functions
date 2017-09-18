#include <fem.hpp>

namespace f2c {
#include <f2c.h>
}

extern "C" /* Subroutine */ int zbesi_(f2c::doublereal *zr, f2c::doublereal *zi, const f2c::doublereal *fnu, const f2c::integer *kode, const f2c::integer *n, f2c::doublereal *cyr, f2c::doublereal *cyi, f2c::integer * nz, f2c::integer *ierr);
extern "C" /* Subroutine */ f2c::doublereal gamln_(f2c::real *z__, f2c::integer *ierr);

namespace ngfem {

    class SpecialCoefficientFunction_ZBESI : public CoefficientFunction
    {
        shared_ptr<CoefficientFunction> arg;

        f2c::doublereal fnu;
        f2c::integer kode;
        f2c::integer n;

      public:

        SpecialCoefficientFunction_ZBESI(shared_ptr<CoefficientFunction> arg_, f2c::doublereal fnu_=1, f2c::integer kode_=2, f2c::integer n_ =2)
          : CoefficientFunction(1, true),
          arg(arg_), fnu(fnu_), kode(kode_), n(n_)
        {
        }

        double Evaluate (const BaseMappedIntegrationPoint & ip) const override 
        {
            throw Exception("SpecialCoefficientFunction_ZBESI::Evaluate called");
        }


        virtual void Evaluate(const BaseMappedIntegrationPoint & ip, FlatVector<Complex> result) const override
        {
            Complex z = arg->EvaluateComplex(ip);
            f2c::doublereal zr,zi;
            zr = z.real();
            zi = z.imag();

            f2c::doublereal cyr, cyi;
            f2c::integer nz;
            f2c::integer ierr;

            zbesi_(&zr, &zi, &fnu, &kode, &n, &cyr, &cyi, &nz, &ierr);
            if(nz>0) cout << "Number of underflows: " << nz << endl;
            if(ierr>0) cout << "Error: " << ierr << endl;

            result(0) = Complex(cyr, cyi);
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
            f2c::real z = arg->Evaluate(ip);
            f2c::integer ierr;
            double res =  gamln_(&z, &ierr);
            if(ierr>0) cout << "Error: " << ierr << endl;
            return res;
        }

    };
}

#undef min
#undef max
#undef abs
#include <python_ngstd.hpp>
PYBIND11_MODULE(special_functions, m) {
    m.def("Bessel", [] (shared_ptr<ngfem::CoefficientFunction> arg, double fnu, int kode, int n) -> shared_ptr<ngfem::CoefficientFunction>
          {
          return make_shared<ngfem::SpecialCoefficientFunction_ZBESI>(arg, fnu, kode, n);
          }, py::arg("Z"), py::arg("FNU")=0, py::arg("KODE")=1, py::arg("N")=2, py::doc(R"DOCSTRING_(

Input
  Z      - Complex argument
  FNU    - DOUBLE PRECISION initial order, FNU>=0
  KODE   - A parameter to indicate the scaling option
           KODE=1  returns
                   CY(L)=I(FNU+L-1,Z), L=1,...,N
               =2  returns
                   CY(L)=exp(-abs(X))*I(FNU+L-1,Z), L=1,...,N
                   where X=Re(Z)
  N      - Number of terms in the sequence, N>=1
)DOCSTRING_")
    );
    m.def("Gamma", [] (shared_ptr<ngfem::CoefficientFunction> arg) -> shared_ptr<ngfem::CoefficientFunction>
          {
          return make_shared<ngfem::SpecialCoefficientFunction_DGAMLN>(arg);

          });
}
