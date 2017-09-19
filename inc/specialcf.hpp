#ifndef SPECIALCF_HPP_INCLUDED___
#define SPECIALCF_HPP_INCLUDED___
#include <fem.hpp>
#include <tuple>

namespace f2c {
#include <f2c.h>
}

// undef defines in f2c.h
#undef min
#undef max
#undef abs
#include <python_ngstd.hpp>

namespace ngfem {
    template <typename T> class T_SpecialCoefficientFunction;

    template <typename TRes, typename TFirstArg, typename ... TArgs> 
    class T_SpecialCoefficientFunction<TRes(*)(TFirstArg, TArgs...) > : public CoefficientFunction
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
            func(func_), arg(arg_), parameters(TFirstArg(), args...)
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
                return std::apply(func, arguments);
            }
        }

        virtual void Evaluate(const BaseMappedIntegrationPoint & ip, FlatVector<Complex> result) const override
        {
            if constexpr(is_complex)
            {
                Complex z = arg->EvaluateComplex(ip);
                auto arguments = parameters;
                std::get<0>(arguments) = z;
                result(0) = std::apply(func, arguments);
            }
            
        }

        template<typename ...TPyArgs>
        static void ExportPython( py::module &m, string name, TFunc func, TPyArgs ... py_args )
        {
            m.def(name.c_str(), [func] (shared_ptr<ngfem::CoefficientFunction> arg, TArgs ... args) -> shared_ptr<ngfem::CoefficientFunction>
                  {
                  return make_shared<ngfem::T_SpecialCoefficientFunction<TRes(*)(TFirstArg, TArgs...) >>(arg, func, args...);
                  },  py_args...
            );
        }
    };

}


template<typename TFunc, typename ...TPyArgs>
static void ExportPythonSpecialCF( py::module &m, string name, TFunc func, TPyArgs ... py_args ) {
    ngfem::T_SpecialCoefficientFunction<TFunc>::ExportPython(m, name, func, py_args...);
}
#endif // SPECIALCF_HPP_INCLUDED___
