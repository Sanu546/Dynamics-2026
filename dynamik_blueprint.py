import sympy as sp
import matplotlib.pyplot as plt
import sympy.plotting as plot
from sympy import Inverse
import numpy as np
from sympy import latex

class DynamikBlueprint:
    def __init__(self):
        pass
# function to calculate jacobian matrix 
    def jacobian(self,matrix_koordates,matrix_variables):
        
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

    # function to calculate velocity and acceleration from position vector
    def posTOacc(self,pos):
        t = sp.symbols('t')
        vel = sp.diff(pos,t)
        acc = sp.diff(vel,t)
        return vel, acc


    def rotation_matrix_extrinsic(self,angles, order='none'):
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
            R = self.order_rotation_matrix(i, sp.rad(angles[x]))
            listorder.append(R)
        
        R = listorder[0]@listorder[1]@listorder[2]
        
        if R == None:
            
            print("Rotation matrix could not be computed")    
        return R

    def order_rotation_matrix(self,type_of_rotation,angle):
        if type_of_rotation == 'x' or type_of_rotation == 'X':
            R = sp.Matrix([[1, 0, 0],
                            [0, sp.cos(angle), -sp.sin(angle)],
                            [0, sp.sin(angle), sp.cos(angle)]])
        elif type_of_rotation == 'y' or type_of_rotation == 'Y':
            R  = sp.Matrix([[sp.cos(angle), 0, sp.sin(angle)],
                                [0, 1, 0],
                                [-sp.sin(angle), 0, sp.cos(angle)]])
        elif type_of_rotation == 'z' or type_of_rotation == 'Z':
                R = sp.Matrix([[sp.cos(angle), -sp.sin(angle), 0],
                                [sp.sin(angle), sp.cos(angle), 0],
                                [0, 0, 1]])
        else:
                print("Type of rotation not defined or something went wrong in order_rotation_matrix")
        return R

    def total_force_gobale(self, forces, gravity, R):
        f_total = sp.Matrix([0,0,0])
        for f in forces:
            f_total += f
        f_total = R @ f_total + gravity
        return f_total

    def gobal_all_moment_drone(self,forces,R,arm_length):
        lokal_n = sp.Matrix([0,0,0])
        angle = sp.pi*2/len(forces)
        
        for i,f in enumerate(forces):
            j = sp.Matrix([sp.cos(i*angle), sp.sin(i*angle), 0])*arm_length
            lokal_n += j.cross(f)
        global_n = R @ lokal_n
        return global_n

    def recursive_speed_and_acc(self,L, joint_traj, offset ,t):
        """
        Generate a recursive trajectory based on the input trajectory and link lengths.

        Parameters:
        L : list of float
            List of link lengths.
        traj : list of sympy expressions
            List of trajectory expressions for each joint.

        Returns:
        list of sympy expressions
            The resulting recursive trajectory.
        """

        n = len(L)
        if len(joint_traj) != n:
            raise ValueError("Length of traj must match length of L")
        
        # Build positions, velocities and accelerations for all links
        angular_velocity = []
        angular_acceleration = []
        body_frame_vecolity = []
        body_frame_acceleration = []
        com_velocity = []
        com_acceleration = []
        s_cn_arm = []
        
        z = sp.Matrix([[0], [0], [1]])

        cumulative_angle = []
        for i in range(n):
            if i == 0:
                cumulative_angle.append(joint_traj[0] + offset)
            else:
                cumulative_angle.append(cumulative_angle[i-1] + joint_traj[i])
        
        for i in range(n):
            if (i == 0):
                s_c0 = self.s_vec(cumulative_angle[0], L[0]/2)
                s_cn_arm.append(s_c0)
                angular_velocity.append(sp.diff(joint_traj[0], t))
                angular_acceleration.append(sp.diff(angular_velocity[0], t))
                com_velocity.append((angular_velocity[0]*z).cross(s_c0))
                com_acceleration.append(sp.diff(com_velocity[0], t))
                body_frame_vecolity.append(sp.Matrix([[0], [0], [0]]))
                body_frame_acceleration.append(sp.Matrix([[0], [0], [0]]))
                continue
                
            s_ci = self.s_vec(cumulative_angle[i], L[i]/2)
            s_li = self.s_vec(cumulative_angle[i-1], L[i-1])
            s_cn_arm.append(s_ci)

            angular_velocity.append(angular_velocity[i-1] + sp.diff(joint_traj[i], t))
            angular_acceleration.append(sp.diff(angular_velocity[i], t))

            body_frame_vecolity.append(body_frame_vecolity[i-1] + (angular_velocity[i-1]*z).cross(s_li))
            body_frame_acceleration.append(body_frame_acceleration[i-1] + (angular_acceleration[i-1]*z).cross(s_li) + (angular_velocity[i-1]*z).cross((angular_velocity[i-1]*z).cross(s_li)))

            com_velocity.append(body_frame_vecolity[i] + (angular_velocity[i]*z).cross(s_ci))
            com_acceleration.append(body_frame_acceleration[i] + (angular_acceleration[i]*z).cross(s_ci) + (angular_velocity[i]*z).cross((angular_velocity[i]*z).cross(s_ci)))


        return body_frame_vecolity, body_frame_acceleration, com_velocity, com_acceleration, angular_velocity, angular_acceleration, s_cn_arm
    def moment_of_inertia_rod_com(m, L):
        """Calculate the moment of inertia of a rectangular rod about its center.

        Parameters:
        m : float
            Mass of the rod.
        L : float
            Length of the rod.
        Returns:
        float
            The moment of inertia of the rod about its center.
        """
        I = (m*L**2)/12
        return I

    def droneTask(self,angles,order,list_forces,global_force,m,arm,angular_velocity,I):
        R = self.rotation_matrix_extrinsic(angles, order)
        #print("Rotation Matrix:\n", latex(R.evalf()))

        global_inertia = R @ I @ R.T
        #print("Global Inertia Matrix:\n", latex(global_inertia.evalf()))
        
        total_force = self.total_force_gobale(list_forces, global_force, R)
        #print("Total Global Force:\n", latex(total_force.evalf()))
        

        total_moment = self.gobal_all_moment_drone(list_forces,R,arm)
        #print("Total Global Moment:\n", latex(total_moment.evalf()))

        acceleration = total_force/m
        #print("Acceleration:\n", latex(acceleration.evalf()))
        angular_acceleration = global_inertia.evalf().inv() @(total_moment - angular_velocity.cross(global_inertia @ angular_velocity))
        #print("Angular Acceleration:\n", latex(angular_acceleration.evalf()))    
        return R, global_inertia, total_force, total_moment, acceleration, angular_acceleration
    ###! UPDATE
    def s_vec(self,q, link_length):
            return sp.Matrix([[link_length*sp.cos(q)], [link_length*sp.sin(q)], [0]])
def bprint(self, string,value):
        print(f"{string}: ${value}$")



