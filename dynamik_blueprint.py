import sympy as sp
import matplotlib.pyplot as plt
import sympy.plotting as plot
from sympy import Inverse
import numpy as np
from sympy import latex
t = sp.symbols('t')
theta1, theta2 = sp.Function('theta1')(t), sp.Function('theta2')(t)
s = sp.Function('s')(t)


# link length
l1 = 0.8
# mass of link  
m = 2 

#s = sp.symbols('s')
# define coordinates 2 dof Pr
x = s+4*sp.cos(theta1)
y = 3+4*sp.sin(theta1)

# define matrices
matrix_koordates = sp.Matrix([x, y])
matrix_variables = sp.Matrix([s, theta1])


F = sp.Matrix([[5],[-2]])
# function to calculate jacobian matrix 
def jacobian(matrix_koordates,matrix_variables):
    
    # Check if the number of rows match
    if matrix_variables.shape[0] != matrix_koordates.shape[0]:
        print("The matrix of variables must have the same number of rows as the matrix of coordinates.")
        return
    
    # Initialize the Jacobian matrix with zeros of appropriate size
    jacobian_matrix = sp.Matrix.zeros(matrix_koordates.shape[0], matrix_variables.shape[0])
    
    # Calculate the partial derivatives
    for i in range(matrix_koordates.shape[0]):
        for j in range(matrix_variables.shape[0]):
            jacobian_matrix[i,j] = sp.diff(matrix_koordates[i], matrix_variables[j])
        
    return jacobian_matrix

j = (jacobian(matrix_koordates,matrix_variables))
print("Jacobian : \n", latex(j))
print("Determinant of Jacobian : ", latex(j.det()))

print("done jacobian")
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


# function to calculate velocity and acceleration from position vector
def posTOacc(pos):
    vel = sp.diff(pos,t)
    acc = sp.diff(vel,t)
    return vel, acc

veclocity_center, acceleration_center = posTOacc(center_pos)
veclocity_center = veclocity_center.subs(t,0.25).evalf()
acceleration_center = acceleration_center.subs(t,0.25).evalf()
I = 1/12*m*l1**2
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

print("r : ", r)
print("Tau : ", tau)

force_giv = sp.Matrix([[5],[-2],[0]])
center_pos = sp.Matrix([sp.cos(theta), sp.sin(theta),0]) * l1/2
end_pos = sp.Matrix([sp.cos(theta), sp.sin(theta),0]) * l1
vel, acc = posTOacc(center_pos)
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



def rotation_matrix_extrinsic(angles, order='none'):
    """_summary_
    Returns the rotation matrix for given angles and order of rotation
    order define noncapital x, y, z 
    Args:
        angles (_type_): _description_
        order (str, optional): _description_. Defaults to 'xyz'.

    Returns:
        _type_: _description_
    """
    if order == 'none':
        print("No order of rotation defined")
        return None
    alpha = sp.rad(angles[0])
    beta = sp.rad(angles[1])
    gamma = sp.rad(angles[2])
    listorder = []
    for x, i in enumerate(order):
        R = order_rotation_matrix(i, angles[x])
        listorder.append(R)
    
    R = listorder[0]@listorder[1]@listorder[2]
    
    if R == None:
        
        print("Rotation matrix could not be computed")    
    return R

def order_rotation_matrix(type_of_rotation,angle):
    if type_of_rotation == 'x':
     R = sp.Matrix([[1, 0, 0],
                    [0, sp.cos(angle), -sp.sin(angle)],
                    [0, sp.sin(angle), sp.cos(angle)]])
    elif type_of_rotation == 'y':
       R  = sp.Matrix([[sp.cos(angle), 0, sp.sin(angle)],
                        [0, 1, 0],
                        [-sp.sin(angle), 0, sp.cos(angle)]])
    elif type_of_rotation == 'z':
        R = sp.Matrix([[sp.cos(angle), -sp.sin(angle), 0],
                        [sp.sin(angle), sp.cos(angle), 0],
                        [0, 0, 1]])
    else:
        print("Type of rotation not defined or something went wrong in order_rotation_matrix")
    return R

def total_force_gobale(forces, gravity, R):
    f_total = sp.Matrix([0,0,0])
    for f in forces:
        f_total += f
    f_total = R @ f_total + gravity
    return f_total

def gobal_all_moment_drone(forces,R,arm_length):
    lokal_n = sp.Matrix([0,0,0])
    for i,f in enumerate(forces):
        j = sp.Matrix([sp.cos(i*sp.pi/2), sp.sin(i*sp.pi/2), 0])*arm_length
        lokal_n += j.cross(f)
    global_n = R @ lokal_n
    return global_n
    


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

print("drone 5:")
R = rotation_matrix_extrinsic(angles, order='xyz')
print("Rotation Matrix:\n", latex(R.evalf()))
global_inertia = R @ I @ R.T
print("Global Inertia Matrix:\n", latex(global_inertia.evalf()))

gravity = g*m
total_force = total_force_gobale(list_of_forces, gravity, R)
print("Total Global Force:\n", latex(total_force.evalf()))

total_moment = gobal_all_moment_drone(list_of_forces,R,0.4)
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



