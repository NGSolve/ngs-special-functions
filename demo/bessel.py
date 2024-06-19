import ngsolve_special_functions
from ngsolve import *
from netgen.geom2d import unit_square
import scipy.special as sp

ngsglobals.msg_level = 0

mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))


n = 0
while n>=0:
    bessel = ngsolve_special_functions.jv(z=(x-0.5)+1j*(y-0.5), v=n)
    Draw(bessel, mesh, 'bessel', sd=5)
    n = float(input('Set order (<0 to abort): '))

