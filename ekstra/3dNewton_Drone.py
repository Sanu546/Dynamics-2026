import numpy as np
import matplotlib.pyplot as plt

theta1 = np.pi / 4
theta2 = 0.1
theta3 = 0.2

f1 = np.array([0, 0, 3.2])
f2 = np.array([0, 0, 1.4])
f3 = np.array([0, 0, 2.5])
f4 = np.array([0, 0, 3.6])
g = np.array([0, 0, -9.81])
omega = np.array([0.5, 0, 0.1])
m = 5
d = 0.5

I = np.array([[5, 0, 0],
              [0, 5, 0],
                [0, 0, 10]])  
def rotz(theta):
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, -s, 0],
                     [s,  c, 0],
                     [0,  0, 1]])

def roty(theta):
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[ c, 0, s],
                     [ 0, 1, 0],
                     [-s, 0, c]])
def rotx(theta):
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[1, 0,  0],
                     [0, c, -s],
                     [0, s,  c]])

R = rotz(theta1) @ rotx(theta2) @ rotz(theta3)

f_total = 0

for f in [f1, f2, f3, f4]:
    f_total += f

 
f_total = R @ f_total

f_total += g * m 
print("Total Local Force:\n", f_total)  


lokal_n = 0
for i,f in enumerate([f1, f2, f3, f4]):
    print("Local Force:\n", f)
    print(i)
    j = np.array([np.cos(i*np.pi/2), np.sin(i*np.pi/2), 0])*d
    lokal_n += np.cross(j,f)

print("Local Moment:\n", lokal_n)

global_n = R @ lokal_n

a = np.array([[1,2],
                [8,1]])

k = np.array([0, -1])

print("Æ Test:\n", a@k)
print("Local Moment:\n", R)
print("Global Force:\n", f_total)
print("Global Moment:\n", global_n)

acceraltion = f_total / m

# golbal interia matrix
I_global = R @ I @ R.T

angular_acc = np.linalg.inv(I_global) @ (global_n - np.cross(omega, I_global @ omega))

print("Angular Acceleration:\n", angular_acc)
print("Linear Acceleration:\n", acceraltion)
print("Inertia Matrix Global:\n", I_global)
