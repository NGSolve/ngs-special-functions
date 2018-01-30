from ngsolve import *
from netgen.geom2d import unit_square

ngsglobals.msg_level = 0

mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))

import ngsolve.special_functions
gamma = ngsolve.special_functions.gammaln(x*y+2.1)

Draw(gamma, mesh, 'gammaln')
