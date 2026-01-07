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

# Opret en instans af klassen
db = DynamikBlueprint()

t = sp.symbols('t')
theta1, theta2 = sp.Function('theta1')(t), sp.Function('theta2')(t)
s = sp.Function('s')(t)

#s = sp.symbols('s')
# define coordinates 2 dof Pr
x = s+4*sp.cos(theta1)
y = 3+4*sp.sin(theta1)

# define matrices
matrix_koordates = sp.Matrix([x, y])
matrix_variables = sp.Matrix([s, theta1])





F = sp.Matrix([[5],[-2]])
j = (db.jacobian(matrix_koordates,matrix_variables))
print("problem 1:")
print("Jacobian : \n", latex(j))
print("Determinant of Jacobian : ", latex(j.det()))

print("done jacobian")

# define variables for problem 2
# link length
l1 = 0.8
# mass of link  
m = 2 

theta = 3 * sp.sin(sp.pi*t)
s = 3 * sp.sin(sp.pi*t)

# theat1 = sp.symbols('theta1')

# theta1_sym = sp.symbols('theta1')
# x_plot = 0 + 4*sp.cos(theta1_sym)
# y_plot = 3 + 4*sp.sin(theta1_sym)

# sp.plot(x_plot, (theta1_sym, 0, 2*sp.pi), title='x vs theta1', ylabel='x')
# sp.plot(y_plot, (theta1_sym, 0, 2*sp.pi), title='y vs theta1', ylabel='y')


z = sp.Matrix([0,0,1])




# angular velocity
velcity_s = sp.diff(s,t)

# angular acceleration
acceleration_s = sp.diff(velcity_s,t)

# velocity in z direction
v_z = z*velcity_s

# substitute t = 0.25
angular_velocity = velcity_s.subs(t,0.25).evalf()
angular_acceleration = acceleration_s.subs(t,0.25).evalf()

# acceleration in center of link
#acceleration_center = l1/2*angular_acceleration
veclocity_center = l1/2*angular_velocity

center= l1/2
center_pos = center * sp.Matrix([sp.cos(theta), sp.sin(theta)])



veclocity_center, acceleration_center = db.posTOacc(center_pos)
veclocity_center = veclocity_center.subs(t,0.25).evalf()
acceleration_center = acceleration_center.subs(t,0.25).evalf()
I = 1/12*m*l1**2
print("problem 2:")
print("angular Velocity : ", angular_velocity)
print("angular Acceleration : ", latex(angular_acceleration))
print("acceleration in center : ", latex(acceleration_center))
print("velocity in center : ", latex(veclocity_center))
print("Inertia : ", I)

Force = m*acceleration_center
Moment = I*angular_acceleration
print("Force : ", Force)
print("Moment : ", Moment)


sc1 = sp.Matrix([[sp.cos(theta)], [sp.sin(theta)]])

sc1 = sc1.subs(t,0.25).evalf()


r = sc1*l1

tau = F.multiply_elementwise(r)



force_giv = sp.Matrix([[5],[-2],[0]])
center_pos = sp.Matrix([sp.cos(theta), sp.sin(theta),0]) * l1/2
end_pos = sp.Matrix([sp.cos(theta), sp.sin(theta),0]) * l1
vel, acc = db.posTOacc(center_pos)
c = l1
acc = acc.subs(t,0.25).evalf()
# task 5
angle = s.subs(t,0.25).evalf()

moment_giv = end_pos.subs(t,0.25).evalf().cross(force_giv) 
#print("moment_giv : ", moment_giv)

reaktion_force =  (m*acc).evalf()
#print("Reaktion force : ", reaktion_force.subs(t,0.25).evalf())
reaktion_moment = center_pos.subs(t,0.25).cross(reaktion_force).evalf()
#print("Reaktion moment : ", reaktion_moment) 
gravity_force = sp.Matrix([0,-9.81,0])*m
gravity_moment = center_pos.subs(t,0.25).cross(gravity_force).evalf()
#print("Gravity force : ", gravity_force)
#print("Gravity moment : ", gravity_moment)

#print("Angular Acceleration : ", (I*angular_acceleration))
#print("check", moment_giv[2] + reaktion_moment[2] )
total_torque =  (I*angular_acceleration) - moment_giv[2] + reaktion_moment[2] - gravity_moment[2]
print("Total Torque : ", total_torque)


# print(jacobian(matrix_koordates,matrix_variables))
# print(jacobian(matrix_koordates,matrix_variables).det())
# problem 5 drone 2022
# define varibles
f1 = sp.Matrix([0,0,0]) # sæt nul for at teste uden thrust ved f1
f2 = sp.Matrix([0,0,25])
f3 = sp.Matrix([0,0,10])
f4 = sp.Matrix([0,0,20])
list_of_forces = [f1,f2,f3,f4]
m = 4
g = sp.Matrix([0,0,-9.82])
I = sp.Matrix([[3,0,0],
                [0,3,0],
                [0,0,1]])
# xyz rot
angles = sp.Matrix([sp.rad(30),sp.rad(5),sp.rad(3)])
angle_velocities = sp.Matrix([0.1,0.5,0.1])

print("Problem 4:")
R = db.rotation_matrix_extrinsic(angles, order='xyz')
print("Rotation Matrix:\n", latex(R.evalf()))
global_inertia = R @ I @ R.T
print("Global Inertia Matrix:\n", latex(global_inertia.evalf()))

gravity = g*m
total_force = db.total_force_gobale(list_of_forces, gravity, R)
print("Total Global Force:\n", latex(total_force.evalf()))

total_moment = db.gobal_all_moment_drone(list_of_forces,R,0.4)
print("Total Global Moment:\n", latex(total_moment.evalf()))

acceleration = total_force/m
print("Acceleration:\n", latex(acceleration.evalf()))
angular_acceleration = global_inertia.evalf().inv() @(total_moment - angle_velocities.cross(global_inertia @ angle_velocities))
print("Angular Acceleration:\n", latex(angular_acceleration.evalf()))    
print("Done")   
# golbal inertia matrix

print("Problem 3")
# define varibels 
theta = sp.Function('theta')(t)
k = sp.Function('k')(t)
V = sp.symbols('V')

omega = sp.diff(theta,t)
l2 = sp.symbols('l2')
V_c1 = sp.diff(k,t)


m2 = sp.symbols('m2')
m1 = sp.symbols('m1')
g = sp.symbols('g')
d = sp.symbols('d')

I_2 = m2*l2**2/12
V_c2 = sp.Matrix([V_c1,-d,0])+sp.Matrix([0,0,omega]).cross(sp.Matrix([l2/2*sp.sin(theta), l2/2*-sp.cos(theta),0]))
V_c2 = sp.diff(sp.Matrix([sp.sin(theta)*l2/2+k,-sp.cos(theta)*l2/2-d]),t)
V_c2 = V_c2[:2,0]
#V_c2 = V_c2.subs(sp.diff(k,t),V).subs(sp.diff(theta,t),omega)



T1 = 1/2*m1*V_c1**2
T2 = 1/2*m2*(V_c2.dot(V_c2))+1/2*I_2*omega**2

#potential energy
V_E1 = m1*g*0
V_E2 = -m2*g*(d+l2/2*sp.cos(theta))

# Lagrian 
L = T2 - V_E2 + T1 - V_E1

Total_force = sp.diff(sp.diff(L,sp.diff(k,t)),t) - sp.diff(L,k)
Total_moment = sp.diff(sp.diff(L,omega),t) - sp.diff(L,theta)

print("vc2 :", latex(V_c2))
print("T1 :", latex(T1))
print("T2 :", latex(T2.nsimplify()))
print("V_E1 :", latex(V_E1.simplify()))
print("V_E2 :", latex(V_E2.simplify()))
print("Lagrangian :", latex(L.simplify()))
print("Total force :", latex(Total_force.simplify()))
print("Total moment :", latex(Total_moment.simplify()))
## factor().nsimplify()
