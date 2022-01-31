#include <specialcf.hpp>
#include <Faddeeva.hh>

extern "C" {
int zbesi_(f2c::doublereal *zr, f2c::doublereal *zi, const f2c::doublereal *fnu, const f2c::integer *kode, const f2c::integer *n, f2c::doublereal *cyr, f2c::doublereal *cyi, f2c::integer * nz, f2c::integer *ierr);
int zbesj_(f2c::doublereal *zr, f2c::doublereal *zi, const f2c::doublereal *fnu, const f2c::integer *kode, const f2c::integer *n, f2c::doublereal *cyr, f2c::doublereal *cyi, f2c::integer * nz, f2c::integer *ierr);
int zbesk_(f2c::doublereal *zr, f2c::doublereal *zi, const f2c::doublereal *fnu, const f2c::integer *kode, const f2c::integer *n, f2c::doublereal *cyr, f2c::doublereal *cyi, f2c::integer * nz, f2c::integer *ierr);
f2c::doublereal gamln_(f2c::real *z__, f2c::integer *ierr);
int zbesh_(f2c::doublereal *zr, f2c::doublereal *zi, f2c::doublereal *fnu, 
	   f2c::integer *kode, f2c::integer *m, f2c::integer *n,
	   f2c::doublereal *cyr, f2c::doublereal *cyi,
	   f2c::integer *nz, f2c::integer *ierr);	   
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

template <typename TFunc>
Complex THankel( TFunc func, Complex z, double order, int m_, int kode_ )
{
  // input arguments 
  f2c::doublereal fnu = order;
  f2c::integer kode = kode_;
  f2c::integer n = 1;
  f2c::integer m = m_;
  f2c::doublereal zr = z.real();
  f2c::doublereal zi = z.imag();

  // outputs
  f2c::doublereal hr, hi;
  f2c::integer nz;
  f2c::integer ierr;

  func(&zr, &zi, &fnu, &kode, &m, &n, &hr, &hi, &nz, &ierr);
  if(nz>0) cout << "Number of underflows: " << nz << endl;
  if(ierr>0) cout << "Error: " << ierr << endl;

  return Complex(hr, hi);
}

extern "C" {
Complex iv ( Complex z, double order) { return zbesIJK(zbesi_, z, order, 1 ); }
Complex ive( Complex z, double order) { return zbesIJK(zbesi_, z, order, 2 ); }
Complex jv ( Complex z, double order) { return zbesIJK(zbesj_, z, order, 1 ); }
Complex jve( Complex z, double order) { return zbesIJK(zbesj_, z, order, 2 ); }
Complex kv ( Complex z, double order) { return zbesIJK(zbesk_, z, order, 1 ); }
Complex kve( Complex z, double order) { return zbesIJK(zbesk_, z, order, 2 ); }

Complex hankel1 (Complex z, double v) { return THankel(zbesh_, z, v, 1, 1 ); }
Complex hankel1e(Complex z, double v) { return THankel(zbesh_, z, v, 1, 2 ); }
Complex hankel2 (Complex z, double v) { return THankel(zbesh_, z, v, 2, 1 ); }
Complex hankel2e(Complex z, double v) { return THankel(zbesh_, z, v, 2, 2 ); }

  
double gammaln(double x_)
{
  f2c::real x = x_;
  f2c::integer ierr;
  double res =  gamln_(&x, &ierr);
  if(ierr>0) cout << "Error: " << ierr << endl;
  return res;
}


}

namespace ngfem {
  
  Complex wofz(Complex z) {return Faddeeva::w(std::complex<double>(z));}
  Complex erf(Complex z) {return Faddeeva::erf(std::complex<double>(z));}  
  Complex erfc(Complex z) {return Faddeeva::erfc(std::complex<double>(z));}  
  Complex erfcx(Complex z) {return Faddeeva::erfcx(std::complex<double>(z));}  
  Complex erfi(Complex z) {return Faddeeva::erfi(std::complex<double>(z));}  
  Complex dawsn(Complex z) {return Faddeeva::Dawson(std::complex<double>(z));}  
  
}



PYBIND11_MODULE(special_functions, m) {
    py::doc doc_string("Same as in scipy.special");
    ExportPythonSpecialCF(m, "gammaln", gammaln, doc_string);

    py::doc docu("Same as in scipy.special, except that the order of the arguments is swapped.");
    ExportPythonSpecialCF(m, "iv",  iv,  py::arg("z"), py::arg("v")=0, docu);
    ExportPythonSpecialCF(m, "ive", ive, py::arg("z"), py::arg("v")=0, docu);
    ExportPythonSpecialCF(m, "jv",  jv,  py::arg("z"), py::arg("v")=0, docu);
    ExportPythonSpecialCF(m, "jve", jve, py::arg("z"), py::arg("v")=0, docu);
    ExportPythonSpecialCF(m, "kv",  kv,  py::arg("z"), py::arg("v")=0, docu);
    ExportPythonSpecialCF(m, "kve", kve, py::arg("z"), py::arg("v")=0, docu);

    ExportPythonSpecialCF(m, "hankel1",  hankel1,  py::arg("z"), py::arg("v")=0,
			  docu);
    ExportPythonSpecialCF(m, "hankel1e", hankel1e, py::arg("z"), py::arg("v")=0,
			  docu);
    ExportPythonSpecialCF(m, "hankel2", hankel2, py::arg("z"), py::arg("v")=0,
			  docu);
    ExportPythonSpecialCF(m, "hankel2e", hankel2e, py::arg("z"), py::arg("v")=0,
			  docu);

    ExportPythonSpecialCF(m, "wofz",  ngfem::wofz,  py::arg("z"), doc_string);
    ExportPythonSpecialCF(m, "erf",   ngfem::erf,   py::arg("z"), doc_string);
    ExportPythonSpecialCF(m, "erfc",  ngfem::erfc,  py::arg("z"), doc_string);
    ExportPythonSpecialCF(m, "erfcx", ngfem::erfcx, py::arg("z"), doc_string);
    ExportPythonSpecialCF(m, "erfi",  ngfem::erfi,  py::arg("z"), doc_string);
    ExportPythonSpecialCF(m, "dawsn", ngfem::dawsn, py::arg("z"), doc_string);

}
