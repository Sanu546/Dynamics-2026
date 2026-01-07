import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent))

import sympy as sp
import matplotlib.pyplot as plt
import sympy.plotting as plot
from sympy import Inverse
import numpy as np
from sympy import latex
from dynamik_blueprint import DynamikBlueprint


db = DynamikBlueprint()

# problem 1

#define variables for problem 1
t = sp.symbols('t') 
theta1 = sp.Function('theta1')(t)
theta2 = sp.Function('theta2')(t)

s = sp.Function('s')(t)

x = 5*sp.cos(theta1)+3*sp.cos(theta1+theta2)
y = 5*sp.sin(theta1)+3*sp.sin(theta1+theta2)
z = 0.5*s+3

# jacobian matrix
j = db.jacobian(sp.Matrix([x,y,z]), sp.Matrix([theta1, theta2, s]))

print("problem 2:")
print("Jacobian:\n", latex(j))
print("Determinant of Jacobian : ", latex(j.det()))


# problem 2
# define variables for problem 2
l1 = 0.5
l2 = 0.3
t = sp.symbols('t')
# in task they use 3.141
theta1 = 3 * sp.sin(3.141*t+1.1)
theta2 = 2 * sp.sin(3.141*t+0.4)

OA_angular_velocity = theta1.diff(t)
AB_angular_velocity = (theta1+theta2).diff(t)



print("Angular velocity OA:\n", latex(OA_angular_velocity.subs(t,0.1)))
print("Angular velocity AB:\n", latex(AB_angular_velocity.subs(t,0.1)))