import sympy as sp
import numpy as np
    
    
t = sp.symbols('t')

omega = sp.Matrix([sp.Function('5')(t),0,sp.Function('2')(t)])
omega = sp.Matrix([5,0,3])

Ix = 10
Iy = 10
Iz = 5

I = sp.Matrix([ [Ix,0,0],
              [0,Iy,0],
              [0,0,Iz]])

domega_dt = sp.Matrix([sp.diff(omega[0],t),
                       sp.diff(omega[1],t),
                       sp.diff(omega[2],t)])
croos = omega.cross( I @ omega)
torque = I * domega_dt + croos
sp.pprint(croos)


