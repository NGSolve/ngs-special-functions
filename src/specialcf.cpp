#include <specialcf.hpp>

extern "C" {
int zbesi_(f2c::doublereal *zr, f2c::doublereal *zi, const f2c::doublereal *fnu, const f2c::integer *kode, const f2c::integer *n, f2c::doublereal *cyr, f2c::doublereal *cyi, f2c::integer * nz, f2c::integer *ierr);
int zbesj_(f2c::doublereal *zr, f2c::doublereal *zi, const f2c::doublereal *fnu, const f2c::integer *kode, const f2c::integer *n, f2c::doublereal *cyr, f2c::doublereal *cyi, f2c::integer * nz, f2c::integer *ierr);
int zbesk_(f2c::doublereal *zr, f2c::doublereal *zi, const f2c::doublereal *fnu, const f2c::integer *kode, const f2c::integer *n, f2c::doublereal *cyr, f2c::doublereal *cyi, f2c::integer * nz, f2c::integer *ierr);
f2c::doublereal gamln_(f2c::real *z__, f2c::integer *ierr);
}

template <typename TFunc>
Complex zbesIJK( TFunc func, Complex z, double order, int kode_ )
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

  func(&zr, &zi, &fnu, &kode, &n, &cyr, &cyi, &nz, &ierr);
  if(nz>0) cout << "Number of underflows: " << nz << endl;
  if(ierr>0) cout << "Error: " << ierr << endl;

  return Complex(cyr, cyi);
}

Complex iv ( Complex z, double order) { return zbesIJK(zbesi_, z, order, 1 ); }
Complex ive( Complex z, double order) { return zbesIJK(zbesi_, z, order, 2 ); }
Complex jv ( Complex z, double order) { return zbesIJK(zbesj_, z, order, 1 ); }
Complex jve( Complex z, double order) { return zbesIJK(zbesj_, z, order, 2 ); }
Complex kv ( Complex z, double order) { return zbesIJK(zbesk_, z, order, 1 ); }
Complex kve( Complex z, double order) { return zbesIJK(zbesk_, z, order, 2 ); }

double gamln(double x_)
{
  f2c::real x = x_;
  f2c::integer ierr;
  double res =  gamln_(&x, &ierr);
  if(ierr>0) cout << "Error: " << ierr << endl;
  return res;
}

#undef min
#undef max
#undef abs
#include <python_ngstd.hpp>
PYBIND11_MODULE(special_functions, m) {
    ExportPythonSpecialCF(m, "Gamma", gamln);

    py::doc docu = "Same as scipy.special.{name}, except that the order of the arguments is swapped.";
    ExportPythonSpecialCF(m, "iv",  iv,  py::arg("z"), py::arg("order")=0, docu);
    ExportPythonSpecialCF(m, "ive", ive, py::arg("z"), py::arg("order")=0, docu);
    ExportPythonSpecialCF(m, "jv",  jv,  py::arg("z"), py::arg("order")=0, docu);
    ExportPythonSpecialCF(m, "jve", jve, py::arg("z"), py::arg("order")=0, docu);
    ExportPythonSpecialCF(m, "kv",  kv,  py::arg("z"), py::arg("order")=0, docu);
    ExportPythonSpecialCF(m, "kve", kve, py::arg("z"), py::arg("order")=0, docu);
}
