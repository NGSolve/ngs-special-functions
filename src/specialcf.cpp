#include <specialcf.hpp>
#include <boost/math/special_functions/gamma.hpp>

extern "C" /* Subroutine */ int zbesi_(f2c::doublereal *zr, f2c::doublereal *zi, const f2c::doublereal *fnu, const f2c::integer *kode, const f2c::integer *n, f2c::doublereal *cyr, f2c::doublereal *cyi, f2c::integer * nz, f2c::integer *ierr);

// extern "C" /* Subroutine */ f2c::doublereal gamln_(f2c::real *z__, f2c::integer *ierr);

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

// double gamln(double x_)
// {
//   f2c::real x = x_;
//   f2c::integer ierr;
//   double res =  gamln_(&x, &ierr);
//   if(ierr>0) cout << "Error: " << ierr << endl;
//   return res;
// }

#undef min
#undef max
#undef abs

namespace bm = boost::math;
#include <python_ngstd.hpp>
PYBIND11_MODULE(special_functions, m) {
  ExportPythonSpecialCF(m,"Gamma",bm::tgamma<double>);
  
  ExportPythonSpecialCF(m, "Bessel", zbesi,
              py::arg("z"), py::arg("order")=0, py::arg("kode")=1, py::doc(R"DOCSTRING_(
Complex Bessel functions
Input
  z      - Complex argument
  order  - DOUBLE PRECISION order>=0
  kode   - A parameter to indicate the scaling option
           kode=1  returns
                   CY=I(order,z)
               =2  returns
                   CY=exp(-abs(X))*I(order,z)
                   where X=Re(z)
)DOCSTRING_"));

}
