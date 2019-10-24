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
import numpy as np
from tp1_modificado import *
# In[]
# Definición de parámetros (en mm)

a1 = 70
a2 = 360
d4 = 380

tol = 1e-6


# In[]
# Defino una función que me devuelve la matriz homogenea que representa 
# la rototraslación a partir del criterio Denavit-Hartenberg 

# ángulos deben estar en radianes
def DH_hom_mat(theta,d,a,alpha):
    
    # Se ponen en rango [-180,180)los ángulos
    #theta = get_in_range(theta)
    #alpha = get_in_range(alpha)
    
    
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
    A[absolute(A)<tol] = 0
    return A

# In[]
    
""" 
Debug DH_hom_mat
    theta = pi/2
    alpha = pi
    a =1
    d =1
"""
# In[]
# Problema directo: angulos deben estar en grados

def pos_prob_dir(t1, t2, t3, t4, t5, t6):    
    # Se ponen en rango [-180,180)los ángulos
    #t1 = get_in_range(t1);t2 = get_in_range(t2);t3 = get_in_range(t3);
    #t4 = get_in_range(t4);t5 = get_in_range(t5);t6 = get_in_range(t6);
    
    t1_act =t1
    t4_act = t4
    [t1,t2,t3,t4,t5, t6] = array([t1,t2,t3,t4,t5, t6])*pi/180
    
    A_1_0  = DH_hom_mat (t1,0,a1, -pi/2)
        
    A_2_1  = DH_hom_mat (t2,0,a2, 0)
    
    A_3_2  = DH_hom_mat (t3,0,0, pi/2)
        
    A_4_3  = DH_hom_mat (t4,d4,0, -pi/2)
        
    A_5_4  = DH_hom_mat (t5,0,0, pi/2)
        
    A_6_5  = DH_hom_mat (t6,0,0, 0)   
    
    A = linalg.multi_dot([A_1_0, A_2_1, A_3_2, A_4_3, A_5_4, A_6_5])
    
    g1 = 1 if sign(d4*sin(t2+t3) +a2*cos(t2) + a1) >= 0 else -1
    g2 = 1 if sign(cos(t3)) >= 0 else -1
    g3 = 1 if sign(t5) >= 0 else -1
    
    
    return (A,[g1,g2,g3],[t1_act,t4_act])

# In[]
# problema iRMat2Eulnverso: los ángulos se devuelven en grados

def pos_prob_inv(A,g1,g2,g3,t1_act=0,t4_act=0):

    px = A[0,3]
    py = A[1,3]
    pz = A[2,3]
    
    
    # cálculo de theta 1
    t1_act = t1_act*pi/180
    t1 = arctan2(g1*py,g1*px) if (absolute(px)> tol or absolute(py)> tol ) else t1_act
    
    # cálculo de theta 3
    s3 = ((px*cos(t1) + py*sin(t1) - a1)**2 + pz**2 -d4**2 -a2**2 ) / (2*a2*d4)
    #s3 = np.round( s3 ,cant_de_decimales)
    if absolute(s3) <= 1:
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
    
    # Mt2 y lt2 son las matrices M y l respectivamente del sistema lineal
    Mt2 = array([[d4*cos(t3), d4*sin(t3) + a2 ],[d4*sin(t3) + a2 , -d4*cos(t3) ]])
    lt2 = array([px*cos(t1) + py*sin(t1) -a1, -pz])
    [s2,c2] = linalg.solve(Mt2,lt2) # Esta función resuelve el sistema lineal (no usa cálculo simbólico)
    t2 = arctan2(s2,c2)
    
    
    # Cálculo de theta 4, 5 y 6
    
    # Calculo A_3_0 canlculando las primeras 3 matrices homogeneas y multiplicándolas
    A_1_0  = DH_hom_mat (t1,0,a1, -pi/2)    
    A_2_1  = DH_hom_mat (t2,0,a2, 0)
    A_3_2  = DH_hom_mat (t3,0,0, pi/2)
    
    A_3_0 = linalg.multi_dot([A_1_0, A_2_1, A_3_2])
    A_3_0[absolute(A_3_0)<tol] =0
    # Obtengo R_3_0 y la matriz de rotación de A
    R_3_0 = A_3_0[0:3,0:3]
    R = A[0:3,0:3]
    # Obtengo la matriz a la que es igual R_6_3 ya que estoy multiplicando
    # la inversa de R_3_0 (.T aplica la operación de transponer) y R
    
    R_p = matmul(R_3_0.T,R) # En teoría, R_p = R_6_3 
    
    #print(R_p)
    
    
    t_4_5_6= RMat2Eul(R_p,g3,t4_act) # calculos los ángulos de Euler usango
                                      # g3 que es el signo de t5 y t4_act
                                      # que es el ángulo t4 actual.
    
    # Pongo a t1, t2 y t3 en un vector
    t_1_2_3 = array([t1,t2,t3])*180/pi

    # Junto todos los resultados
    return concatenate([t_1_2_3, t_4_5_6])
    
     

# In[]
"""    
# Pruebas simples
    
#### PRUEBA 1 ####
[A,g,t_act]= pos_prob_dir(0,0,0,0,0,0)
g1 = g[0];g2 = g[1];g3 = g[2];
t1_act=t_act[0];t4_act=t_act[1];
pos_prob_inv(A,*g,*t_act)

#### PRUEBA 2 ####
[A,g,t_act]= pos_prob_dir(0,45,0,0,0,0)
g1 = g[0];g2 = g[1];g3 = g[2];
t1_act=t_act[0];t4_act=t_act[1];
pos_prob_inv(A,*g,*t_act)


#### PRUEBA 3 ####
[A,g,t_act]= pos_prob_dir(0,0,45,0,0,0)
g1 = g[0];g2 = g[1];g3 = g[2];
t1_act=t_act[0];t4_act=t_act[1];
pos_prob_inv(A,*g,*t_act)

#### PRUEBA 4 ####
[A,g,t_act]= pos_prob_dir(0,0,0,45,0,0)
g1 = g[0];g2 = g[1];g3 = g[2];
t1_act=t_act[0];t4_act=t_act[1];
pos_prob_inv(A,*g,*t_act)


#### PRUEBA 5 ####
[A,g,t_act]= pos_prob_dir(0,0,0,0,45,0)
g1 = g[0];g2 = g[1];g3 = g[2];
t1_act=t_act[0];t4_act=t_act[1];
pos_prob_inv(A,*g,*t_act)

#### PRUEBA 6 ####
[A,g,t_act]= pos_prob_dir(0,0,0,0,0,45)
g1 = g[0];g2 = g[1];g3 = g[2];
t1_act=t_act[0];t4_act=t_act[1];
pos_prob_inv(A,*g,*t_act)




### PRUEBA 7 ####       [  0.  45.  90.   0.   0. 135.]
[A,g,t_act]= pos_prob_dir(0,45,90,0,0,135)
g1 = g[0];g2 = g[1];g3 = g[2];
t1_act=t_act[0];t4_act=t_act[1];
pos_prob_inv(A,*g,*t_act)
pos_prob_inv(A,1,-1,-1,*t_act)


### PRUEBA 8 ####       [  15. -180.    0. -135.    0.  -45.]
[A,g,t_act]= pos_prob_dir(15,-180,0,-135,0,-45)
g1 = g[0];g2 = g[1];g3 = g[2];
t1_act=t_act[0];t4_act=t_act[1];
pos_prob_inv(A,*g,*t_act)
[A_2,g_2,t_act_2]= pos_prob_dir(*pos_prob_inv(A,*g,*t_act))
(A -A_2)< tol # Devuelve una matriz con sólo valores True


### PRUEBA 9 ####       [   0.  -45.   45.  135.  180. -180.]
[A,g,t_act]= pos_prob_dir(0,-45,45,135,180,-180)
g1 = g[0];g2 = g[1];g3 = g[2];
t1_act=t_act[0];t4_act=t_act[1];
pos_prob_inv(A,*g,*t_act)




### Prueba TP 1 ###
# Estos ángulos corresponden a los valores que se deben tomar en la figura donde 
# tomaron las ternas

[A,g,t_act]= pos_prob_dir(0,-90,180,0,0,0)
g1 = g[0];g2 = g[1];g3 = g[2];
t1_act=t_act[0];t4_act=t_act[1];
pos_prob_inv(A,*g,*t_act)
pos_prob_inv(A,1,-1,0,*t_act)


### Prueba TP 2 ###
[A,g,t_act]= pos_prob_dir(0,0,0,0,90,0)
g1 = g[0];g2 = g[1];g3 = g[2];
t1_act=t_act[0];t4_act=t_act[1];
pos_prob_inv(A,*g,*t_act)
pos_prob_inv(A,-1,-1,0,*t_act)


### Prueba TP 3 ###
[A,g,t_act]= pos_prob_dir(30,30,30,30,30,30)
g1 = g[0];g2 = g[1];g3 = g[2];
t1_act=t_act[0];t4_act=t_act[1];
pos_prob_inv(A,*g,*t_act)

"""