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
from tp1 import *
# In[]
# Definición de parámetros (en mm)

a1 = 70
a2 = 360
d4 = 380


# In[]
# Defino una función que me devuelve la matriz homogenea que representa 
# la rototraslación a partir del criterio Denavit-Hartenberg 

# ángulos deben estar en radianes
def DH_hom_mat(theta,d,a,alpha):
    R = array(
        [[cos(theta), -sin(theta)*cos(alpha), sin(theta)*sin(alpha)],
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
# Problema directo: angulos deben estar en grados

def pos_prob_dir(t1, t2, t3, t4, t5, t6):    
    
    a1 = 70
    a2 = 360
    d4 = 380

    [t1,t2,t3,t4,t5, t6] = array([t1,t2,t3,t4,t5, t6])*pi/180
    
    A_1_0  = DH_hom_mat (t1,0,a1, -pi/2)
        
    A_2_1  = DH_hom_mat (t2,0,a2, 0)
    
    A_3_2  = DH_hom_mat (t3,0,0, pi/2)
        
    A_4_3  = DH_hom_mat (t4,d4,0, -pi/2)
        
    A_5_4  = DH_hom_mat (t5,0,0, pi/2)
        
    A_6_5  = DH_hom_mat (t6,0,0, 0)   
    
    A = linalg.multi_dot([A_1_0, A_2_1, A_3_2, A_4_3, A_5_4, A_6_5])
    
    g1 = sign(d4*sin(t2+t3) +a2*cos(t2) + a1)
    g2 = sign(cos(t3))
    g3 = sign(t5)
    t4_act = t4*180/pi
    
    return (A,[g1,g2,g3],t4_act)

# In[]
# problema inverso: los ángulos se devuelven en grados

def pos_prob_inv(A,g1,g2,g3,t4_act):
    
    a1 = 70
    a2 = 360
    d4 = 380
    tol = 1e-12
    
    px = A[0,3]
    py = A[1,3]
    pz = A[2,3]
    
    # cálculo de theta 1
    t1 = arctan2(g1*py,g1*px)
    
    # cálculo de theta 3
    s3 = ( (px*cos(t1) + py*sin(t1) - a1)**2 + pz**2 -d4**2 -a2**2 ) / (2*a2*d4)
    if absolute(s3) < 1:
        t3 = arctan2(s3,g2*(1-s3**2)**(0.5)) 
    else:
        print("Posición imposible de obtener!! \n Verificar distancias px, py, pz")
        return [0,0,0,0,0,0]
    
    #cálculo de theta 2
    """
        Sistema lineal de 2 x 2 : Mx = l
        Con M = [d4*cos(t3)         | d4*sin(t3) + a2 ]
                [d4*sin(t3) + a2    | -d4*cos(t3)     ]
                
            x = [ sin(t2)]
                [ cos(t2)]
                
            l = [px*cos(t1) +py*sin(t1) -a2]
                [-pz]
    """
     
    Mt2 = array([[d4*cos(t3), d4*sin(t3) + a2 ],[d4*sin(t3) + a2 , -d4*cos(t3) ]])
    lt2 = array([px*cos(t1) + py*sin(t1) -a1, -pz])
    [s2,c2] = linalg.solve(Mt2,lt2)
    t2 = arctan2(s2,c2)
    
    
    # Cálculo de theta 4, 5 y 6
    
    A_1_0  = DH_hom_mat (t1,0,a1, -pi/2)    
    A_2_1  = DH_hom_mat (t2,0,a2, 0)
    A_3_2  = DH_hom_mat (t3,0,0, pi/2)
    
    A_3_0 = linalg.multi_dot([A_1_0, A_2_1, A_3_2])
    R_3_0 = A_3_0[0:3,0:3]
    R = A[0:3,0:3]
    R_p = matmul(R_3_0.T,R) # En teoría, R_p = R_6_3 
    
    
    """
    ATENCIÓN: tengo que agregar la siguiente linea porque tenía el siguiente
    problema. Cuando probaba la funcion directa con todos los valores de t_i = 0
    excepto t2 (por ejemplo t2 = 45) obtenia un resultado y cuando
    realizaba el paso inverso. t4 y t6 resultaban iguales en módulo pero con 
    signo contrario. Probando encontré que cuando R_p era muy similar a la identidad
    esto es, los elementos de la diagonal eran uno pero los elementos fuera de
    la diagonal eran del orden de 1e-15, este error sucedía. Utilizando los mismos
    valores encontré que si reemplazaba R_p por exactamente identidad, obtenía resultados
    coherentes. Por lo tanto decidí "mandar" a cero cualquier elemento de R_p
    que sea menor a una tolerancia.
    """    
    R_p[absolute(R_p)<tol] = 0
    
    
    t_4_5_6 = RMat2Eul(R_p,g3,t4_act)
    
    t_1_2_3 = array([t1,t2,t3]) * 180/pi
    t_1_2_3[absolute(t_1_2_3)<tol]  = 0 # agrego esto, no porque exista un error
                                        # sino para que sea más fácil visualizar
                                        # el resultado.
    
    theta_vec = concatenate([t_1_2_3, t_4_5_6])
    
    return theta_vec

# In[]
    
# Pruebas simples
    
[A,g,t4_act] = pos_prob_dir(0,0,0,0,0,0)
pos_prob_inv(A,*g,t4_act)


[A,g,t4_act]= pos_prob_dir(-135,0,0,0,0,0)
pos_prob_inv(A,*g,t4_act)



t1 = 0; t2 = 45; t3 = 0; t4 = 0; t5 = 0; t6 = 0;
[A,g,t4_act]= pos_prob_dir(0,45,0,0,0,0)
pos_prob_inv(A,*g,t4_act)


[A,g,t4_act]= pos_prob_dir(0,0,45,45,0,0)
pos_prob_inv(A,*g,t4_act)


[A,g,t4_act]= pos_prob_dir(0,1,-1,45,2,-2)
pos_prob_inv(A,*g,t4_act)

[A,g,t4_act] = pos_prob_dir(0,0,0,0,45,0)
pos_prob_inv(A,*g,t4_act)


[A,g,t4_act]= pos_prob_dir(0,0,0,0,0,45)
pos_prob_inv(A,*g,t4_act)





[A,g,t4_act] = pos_prob_dir(-45,45,45,0,0,0)
pos_prob_inv(A,*g)

[A,g,t4_act] = pos_prob_dir(0,-45,0,0,0,0)
pos_prob_inv(A,*g)

[A,g,t4_act]= pos_prob_dir(45,-45,-45,0,0,0)
pos_prob_inv(A,*g)

