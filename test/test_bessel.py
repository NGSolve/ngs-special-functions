import ngsolve as ngs
from netgen.geom2d import unit_square
import scipy.special as sp
import ngsolve.special_functions as ng_sp
import numpy

ngs.ngsglobals.msg_level = 0
ngs.SetNumThreads(4)

fs_ng = [ng_sp.iv, ng_sp.ive, ng_sp.jv, ng_sp.jve, ng_sp.kv, ng_sp.kve,
         ng_sp.hankel1, ng_sp.hankel1e, ng_sp.hankel2, ng_sp.hankel2e]
fs_py = [sp.iv, sp.ive, sp.jv, sp.jve, sp.kv, sp.kve,
         sp.hankel1, sp.hankel1e, sp.hankel2, sp.hankel2e]


def test_bessel():
    mesh = ngs.Mesh(unit_square.GenerateMesh(maxh=0.2))
    for order in range(10):
        functions = [(f1(ngs.x + 1j*ngs.y, order), f2)
                     for f1, f2 in zip(fs_ng, fs_py)]

        for f_ng, f_py in functions:
            for x in numpy.linspace(0.1, 0.9, 20):
                for y in numpy.linspace(0.1, 0.9, 20):
                    z = x+1j*y
                    v1 = f_ng(mesh(x, y))
                    v2 = f_py(order, z)
                    error = abs(v1 - v2)
                    if abs(v2) > 1:
                        error = abs(error/v2)
                    assert error < 1e-13


def test_parallel():
    mesh = ngs.Mesh(unit_square.GenerateMesh(maxh=0.2))
    fes = ngs.L2(mesh, order=5, complex=True)
    g1 = ngs.GridFunction(fes)
    g2 = ngs.GridFunction(fes)

    for order in range(10):
        functions = [(f(ngs.x + 1j*ngs.y, order)) for f in fs_ng]

        for f in functions:
            g1.Set(f)
            with ngs.TaskManager():
                g2.Set(f)

            error = ngs.Integrate(g1-g2, mesh)
            assert error == 0j


def test_code_generation():
    mesh = ngs.Mesh(unit_square.GenerateMesh(maxh=0.2))
    for order in range(10):
        functions = [f(ngs.x + 1j*ngs.y, order) for f in fs_ng]

        for f in functions:
            f_err = f-f.Compile(True, wait=True)
            f_err = f_err*ngs.Conj(f_err)
            error = ngs.Integrate(f_err, mesh)
            assert abs(error) < 1e-13


if __name__ == "__main__":
    test_bessel()
    test_parallel()
    test_code_generation()
