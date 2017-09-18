#include <fem.hpp>
#include <tuple>
#include <experimental/tuple>

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

template<typename ...Args> struct FuncTraits{};

template<typename TRes, typename TFirstArg, typename ...TArgs>
struct FuncTraits<TRes(TFirstArg, TArgs...)>
{
  typedef TRes res_type;
  typedef TFirstArg first_arg_type;
  typedef std::tuple<TArgs...> args_tuple_type;

};

namespace ngfem {
    template<typename TRes, typename TFirstArg, typename ...TArgs>
    class T_SpecialCoefficientFunction : public CoefficientFunction
    {
        static constexpr bool is_complex = std::is_same<Complex, TRes>::value;
        static constexpr bool is_double = std::is_same<double, TRes>::value;
        static_assert(is_complex || is_double, "Result has to be either Complex or double");

        typedef TRes (*TFunc)(TFirstArg, TArgs...);

        TFunc func;
        shared_ptr<CoefficientFunction> arg;
        std::tuple<TFirstArg, TArgs...> parameters;

      public:
        T_SpecialCoefficientFunction( shared_ptr<CoefficientFunction> arg_, TFunc func_, TArgs ... args )
          : CoefficientFunction(1, is_complex),
            func(func_), arg(arg_), parameters({}, args...)
        {
        }

        double Evaluate (const BaseMappedIntegrationPoint & ip) const override 
        {
            if constexpr(is_complex)
              throw Exception("T_SpecialCoefficientFunction::Evaluate called for complex function");
            else {
                double z = arg->Evaluate(ip);
                auto arguments = parameters;
                std::get<0>(arguments) = z;
                return std::experimental::apply(func, arguments);
            }
        }

        virtual void Evaluate(const BaseMappedIntegrationPoint & ip, FlatVector<Complex> result) const override
        {
            if constexpr(is_complex)
            {
                Complex z = arg->EvaluateComplex(ip);
                auto arguments = parameters;
                std::get<0>(arguments) = z;
                result(0) = std::experimental::apply(func, arguments);
            }
            
        }

    };

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
    m.def("Gamma", [] (shared_ptr<ngfem::CoefficientFunction> arg) -> shared_ptr<ngfem::CoefficientFunction>
          {
          return make_shared<ngfem::T_SpecialCoefficientFunction<double, double>>(arg, gamln);
          });

    m.def("Bessel", [] (shared_ptr<ngfem::CoefficientFunction> arg, double order, int kode) -> shared_ptr<ngfem::CoefficientFunction>
          {
          return make_shared<ngfem::T_SpecialCoefficientFunction<Complex, Complex, double, int>>(arg, zbesi, order, kode);
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
}
