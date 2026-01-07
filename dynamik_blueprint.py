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
            R = self.order_rotation_matrix(i, angles[x])
            listorder.append(R)
        
        R = listorder[0]@listorder[1]@listorder[2]
        
        if R == None:
            
            print("Rotation matrix could not be computed")    
        return R

    def order_rotation_matrix(self,type_of_rotation,angle):
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

    def total_force_gobale(self, forces, gravity, R):
        f_total = sp.Matrix([0,0,0])
        for f in forces:
            f_total += f
        f_total = R @ f_total + gravity
        return f_total

    def gobal_all_moment_drone(self,forces,R,arm_length):
        lokal_n = sp.Matrix([0,0,0])
        for i,f in enumerate(forces):
            j = sp.Matrix([sp.cos(i*sp.pi/2), sp.sin(i*sp.pi/2), 0])*arm_length
            lokal_n += j.cross(f)
        global_n = R @ lokal_n
        return global_n




