from ngsolve import *
from netgen.geom2d import unit_square

ngsglobals.msg_level = 0

mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))

import ngsolve.special_functions

arg = x+0.1+1j*(y+0.1) 
bessel = ngsolve.special_functions.Bessel(Z=0.00001*arg, FNU=0.1, KODE=1, N=1)

Draw(bessel, mesh, 'bessel')

