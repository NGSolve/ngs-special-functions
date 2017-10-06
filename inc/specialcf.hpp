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

template<typename T>
inline string GetTypeName()
{
    if constexpr(std::is_same_v<T, double>) return "double";
    if constexpr(std::is_same_v<T, Complex>) return "Complex";
}

namespace ngfem {

    template <typename T, typename TFunc>
    void IterateTuple(const T & t, const TFunc &f)
    {
        constexpr size_t n = std::tuple_size_v<T>;
        static_assert(n < 10, "Increase number of lines in IterateTuple.");
        if constexpr(n>0) f(std::get<0>(t));
        if constexpr(n>1) f(std::get<1>(t));
        if constexpr(n>2) f(std::get<2>(t));
        if constexpr(n>3) f(std::get<3>(t));
        if constexpr(n>4) f(std::get<4>(t));
        if constexpr(n>5) f(std::get<5>(t));
        if constexpr(n>6) f(std::get<6>(t));
        if constexpr(n>7) f(std::get<7>(t));
        if constexpr(n>8) f(std::get<8>(t));
        if constexpr(n>9) f(std::get<9>(t));
    }


    template <typename T> class T_SpecialCoefficientFunction;

    template <typename TRes, typename TFirstArg, typename ... TArgs> 
    class T_SpecialCoefficientFunction<TRes(*)(TFirstArg, TArgs...) > : public CoefficientFunction
    {
        static constexpr int nargs = sizeof...(TArgs)+1;
        static constexpr bool is_complex = std::is_same<Complex, TRes>::value;
        static constexpr bool is_double = std::is_same<double, TRes>::value;
        static_assert(is_complex || is_double, "Result has to be either Complex or double");

        typedef TRes (*TFunc)(TFirstArg, TArgs...);

        TFunc func;
        string name;
        shared_ptr<CoefficientFunction> arg;
        std::tuple<TFirstArg, TArgs...> parameters;

      public:
      T_SpecialCoefficientFunction( shared_ptr<CoefficientFunction> arg_, TFunc func_, string name_, TArgs ... args )
          : CoefficientFunction(1, is_complex),
            func(func_), name(name_), arg(arg_), parameters(TFirstArg(), args...)
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

      virtual void TraverseTree (const function<void(CoefficientFunction&)> & func) override
      {
        arg->TraverseTree (func);
        func(*this);
      }

      virtual Array<CoefficientFunction*> InputCoefficientFunctions() const override
      { return Array<CoefficientFunction*>({ arg.get() }); }

      virtual void GenerateCode(Code & code, FlatArray<int> inputs, int index) const override
      {
        string params = "";
        int ii = 0;
        IterateTuple( parameters, [&](auto p) {
                     if(ii++) // skip first argument
                     params += ", " + ToString(p);
               });


        if(code.is_simd==false && code.deriv==0)
        {
            string declaration = string("extern ") + GetTypeName<TRes>() + " " + name + " (";
            int ii = 0;
            IterateTuple( parameters, [&](auto p) {
                     declaration += (ii++ ? ", " : "") + GetTypeName<decltype(p)>();
            });
            declaration += ");\n";
            if(code.top.find(declaration) == string::npos)
                code.top += declaration;
            code.AddLinkFlag(SPECIALCF_LIBRARY_NAME);
        }

        if(code.is_simd)
        {
            auto tmp_real = Var("tmp_real",index);
            auto tmp_imag = Var("tmp_imag",index);
            code.body += tmp_real.Declare( "SIMD<double>" );
            if(is_complex) code.body += Var("tmp_imag",index).Declare( "SIMD<double>" );

            for (int i : Range(SIMD<double>::Size()))
            {
              auto tmp = Var("tmp",index,i);
              string comp = "["+ToLiteral(i)+"]";
              if(is_complex)
              {
                  auto var = Var(inputs[0]).S();
                  code.body += tmp.Assign( name+"(Complex("+var+".real()"+comp+","+var+".imag()"+comp+")"+params+")" );
                  code.body += tmp_real.S() + comp + " = " + tmp.S() + ".real();\n";
                  code.body += tmp_imag.S() + comp + " = " + tmp.S() + ".imag();\n";
              }
              else
              {
                  code.body += tmp.Assign( name+"("+Var(inputs[0]).S()+comp+params+")" );
                  code.body += tmp_real.S() + comp + " = " + tmp.S() + ";\n";
              }

            }

            if(is_complex)
                code.body += Var(index).Assign( "SIMD<Complex>("+tmp_real.S() + "," + tmp_imag.S() + ")" );
            else
                code.body += Var(index).Assign( tmp_real.S() );
        }
        else
            code.body += Var(index).Assign( name+"("+Var(inputs[0]).S()+params+")" );

      }

        template<typename ...TPyArgs>
        static void ExportPython( py::module &m, string name, TFunc func, TPyArgs ... py_args )
        {
          m.def(name.c_str(), [func,name] (shared_ptr<ngfem::CoefficientFunction> arg, TArgs ... args) -> shared_ptr<ngfem::CoefficientFunction>
                  {
                    return make_shared<ngfem::T_SpecialCoefficientFunction<TRes(*)(TFirstArg, TArgs...) >>(arg, func, name, args...);
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
