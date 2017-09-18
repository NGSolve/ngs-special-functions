from ngsolve import *
from netgen.geom2d import unit_square

ngsglobals.msg_level = 0

mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))

import ngsolve.special_functions

n = 0
while n>=0:
    bessel = ngsolve.special_functions.Bessel(z=(x-0.5)+1j*(y-0.5), order=n, kode=1)
    Draw(bessel, mesh, 'bessel', sd=5)
    n = float(input('Set order (<0 to abort): '))

