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
L = [l1, l2]
t = sp.symbols('t')
# in task they use 3.141
theta1 = 3 * sp.sin(3.141*t+1.1)
theta2 = 2 * sp.sin(3.141*t+0.4)

OA_angular_velocity = theta1.diff(t)
AB_angular_velocity = (theta1+theta2).diff(t)


t_def = 0.1
theta_list = [theta1, theta2]
bv_list, ba_list, com_vel_list, com_acc_list, av_list, aa_list, s_cn_list = db.recursive_speed_and_acc(L,theta_list,-1/2*sp.pi ,t)
print("Angular velocity OA:\n", latex(OA_angular_velocity.subs(t,0.1)))
print("Angular velocity AB:\n", latex(AB_angular_velocity.subs(t,0.1)))

print("Anular velocity OA:\n", av_list[0].subs(t,t_def).evalf())
print("Anular velocity AB:\n", av_list[1].subs(t,t_def).evalf())

print("Anular acceleration OA:\n", aa_list[0].subs(t,t_def).evalf())
print("Anular acceleration AB:\n", aa_list[1].subs(t,t_def).evalf())

print("COM velocity C1:\n", latex(com_vel_list[0].subs(t,t_def).evalf()))
print("COM velocity C2:\n", latex(com_vel_list[1].subs(t,t_def).evalf()))
    
print("COM acceleration C1:\n", latex(com_acc_list[0].subs(t,t_def).evalf()))
print("COM acceleration C2:\n", latex(com_acc_list[1].subs(t,t_def).evalf()))

F_A =sp.Matrix([[sp.symbols('F_Ax')], [sp.symbols('F_Ay')]])
F_O =sp.Matrix([[sp.symbols('F_Ox')], [sp.symbols('F_Oy')]])
m1 = 3
m2 = 2
I1 = 0.5
I2 = 0.3

F_A = m2*com_acc_list[1]+m2*sp.Matrix([[0],[9.81],[0]])
F_O = m1*com_acc_list[0]+m1*sp.Matrix([[0],[9.81],[0]])+F_A
sc1 = sp.Matrix([[ l1*sp.cos(theta1-sp.pi/2)/2],[ l1*sp.sin(theta1-sp.pi/2)/2],[0]])
sc_1 = sp.Matrix([[ l1*sp.cos(theta1-sp.pi/2)],[ l1*sp.sin(theta1-sp.pi/2)],[0]])
sc2 = sc_1+sp.Matrix([[l2*sp.cos(theta1+theta2-sp.pi/2)/2],[l2*sp.sin(theta1+theta2-sp.pi/2)/2],[0]])
print("I2 ",I2*aa_list[1].subs(t,t_def).evalf())
print("sc2 cross F_A ",(s_cn_list[1].cross(F_A)).subs(t,t_def).evalf())
T_2 = I2*aa_list[1]+(s_cn_list[1].cross(F_A))[2]
T_1 = I2*aa_list[0]+(s_cn_list[0].cross(F_O))[2]+(s_cn_list[1].cross(F_A))[2]+T_2

print("Force at A:\n", latex(F_A.subs(t,t_def).evalf()))
print("Force at O:\n", latex(F_O.subs(t,t_def).evalf()))
print("Torque at 2:\n", latex(T_2.subs(t,t_def).evalf()))
print("Torque at 1:\n", latex(T_1.subs(t,t_def).evalf()))

# problem 3
print("problem 3:")
t = sp.symbols('t') 
theta,h = sp.Function('theta')(t),sp.Function('h')(t)
l,m1,m2 = sp.symbols('l'),sp.symbols("m1"),sp.symbols("m2")
theta1 = sp.Function('theta1')(t)
I = m2*l2**2/12
sc1 = sp.Matrix([[0],[h],[0]])
sc2 = sc1 + l/2*sp.Matrix([[sp.cos(theta)],[sp.sin(theta)],[0]])

v_c2 = sp.diff(sc2,t)
print(sc2)
print("VC2",latex(v_c2))






# problem 4

print("problem 4:")
m = 6
print("Mass :", )
I = sp.Matrix([[4,0,0],
                [0,4,0],
                [0,0,2]])
angels = [45, 0, 5]
omega = sp.Matrix([[0.1],[0.1],[0.7]])
F_1 = sp.Matrix([[0],[0],[10]])
F_2 = sp.Matrix([[0],[0],[25]])
F_3 = sp.Matrix([[0],[0],[5]])
F_4 = sp.Matrix([[0],[0],[30]])
list_of_forces = [F_1, F_2, F_3, F_4]
arm = 0.4
p = sp.Matrix([[6],[-3],[-30]])
g = sp.Matrix([[0],[0],[-9.81]])
gravyti_force = g*m
global_force = gravyti_force+p


R, global_inertia, total_force, total_moment, acceleration, angular_acceleration = db.droneTask(angels,"zyx",list_of_forces,global_force,m,arm,omega,I)
print("Roation Matrix:\n", latex(R.evalf()))
print("Global Inertia Matrix:\n", latex(global_inertia.evalf()))
print("acceleration:\n", latex(acceleration.evalf()))
print("angular_acceleration:\n", latex(angular_acceleration.evalf()))