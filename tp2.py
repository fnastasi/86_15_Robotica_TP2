#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 12:31:14 2019

@author: franco
"""

"""
El carácter * aplicado a un vector "desarma el vector". 
Solo se puede aplicar al parámetro de una función. 
Ejemplo: x = [1,2,3] => fn(*x) == fn(1,2,3)
"""

"""
Problema directo 

Ejemplo de uso:
    theta_vec = [45,45,45,30,30,30]
    R, G = PosDir(*theta_vec)
"""
# In[]

from numpy import *

# In[]
# Definición de parámetros (en mm)

a1 = 70
a2 = 360
d4 = 380


# In[]
# Defino una función que me devuelve la matriz homogenea que representa 
# la rototraslación a partir del criterio Denavit-Hartenberg 

def DH_hom_mat(theta,d,a,alpha):
    R = array(
        [[cos(theta), -sin(theta)*cos(theta), sin(theta)*sin(alpha)],
        [sin(theta), cos(theta)*cos(alpha), -cos(theta)*sin(alpha)],
        [0,sin(alpha) , cos(alpha)]])
    P = array([[a*cos(theta),
               a*sin(theta),
               d]]).T
    A = concatenate((R,P),axis = 1) # Junta R y P -> [R | P]
    
    # Junta Ry P con el vector [0001] -> [  R    | P ]
    #                                    [ 0 0 0 | 1 ]
    
    A = concatenate((A,array([[0,0,0,1]])), axis = 0) 
    return A

# In[]
    
# Debug DH_hom_mat
    theta = 0
    alpha = 0
    a =1
    d =1

# In[]
        
def pos_prob_dir(t1, t2, t3, t4, t5, t6):    
    
    a1 = 70
    a2 = 360
    d4 = 380


    A_1_0  = DH_hom_mat (t1,0,a1, -pi/2)
        
    A_2_1  = DH_hom_mat (t2/2,0,a2, 0)
    
    A_3_2  = DH_hom_mat (t3,0,0, pi/2)
        
    A_4_3  = DH_hom_mat (t4,d4,0, -pi/2)
        
    A_5_4  = DH_hom_mat (t5,0,0, pi/2)
        
    A_6_5  = DH_hom_mat (t6,0,0, 0)   
    
    A = linalg.multi_dot([A_1_0, A_2_1, A_3_2, A_4_3, A_5_4, A_6_5])
    
    g1 = sign(d4*sin(t2+t3) +a2*cos(t2) + a1)
    g2 = sign(cos(t3))
    g3 = 0
    
    return (A,[g1,g2,g3])

# In[]
def pos_prob_inv(A,g1,g2,g3):
    
    
    return theta_vec
