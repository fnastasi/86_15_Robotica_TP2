#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 18:16:41 2019

@author: franco
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 08:39:11 2019

@author: franco
"""

# In[]
# Se importa la biblioteca numpy para aplicar operaciones similares a las que se usan en matlab
from numpy import *
import numpy as np


# In[]

tol = 1e-6

# In[]

"""
Eul2RMat es la función que resuelve el problema directo de los ángulos de Euler
Los parámetros de entrada son los ángulos tita, phi y psi
Ejemplo de como utilizar la función para phi = 45°, tita = 45° y psi =0° :
    (R,sg_tita) = Eul2RMat(phi =45,tita=45,psi=0 )
    (R,sg_tita) = Eul2RMat(45,45,0 )

Precondición:
    - Los ángulos deben ingresarse en grado
    - phi y psi pueden encontrarse en el rango [-pi, pi]
    - tita debe encontrarse en el rango [0, pi]
"""


def Eul2RMat(phi=0, tita=0,psi=0):
    # los ángulos se pasan a radianes
    phi = phi*pi/180
    tita = tita*pi/180
    psi = psi*pi/180
    
    # Los ángulos se ponen en rango [-180,180)
    #phi = get_in_range(phi)*pi/180
    #tita = get_in_range(tita)*pi/180
    #psi = get_in_range(psi)*pi/180
    
    # Se escriben las matrices de rotación y con la función array se convierten en elementos
    # de numpy
    R_z_phi = [[cos(phi), -sin(phi), 0],[sin(phi), cos(phi), 0],[0,0,1]]
    R_z_phi = array(R_z_phi)
    #R_z_phi = np.round(array(R_z_phi),cant_de_decimales)
    #print(R_z_phi)
    
    R_y_tita = [[cos(tita), 0, sin(tita)], [0,1,0] , [-sin(tita), 0, cos(tita)]]
    R_y_tita = array(R_y_tita)
    #R_y_tita = np.round(array(R_y_tita),cant_de_decimales)
    #print(R_y_tita)
    
    R_z_psi = [[cos(psi), -sin(psi), 0],[sin(psi), cos(psi), 0],[0,0,1]]
    R_z_psi =array(R_z_psi)
    #R_z_psi = np.round(array(R_z_psi),cant_de_decimales)
    
    #print(R_z_psi)
    
    # Matriz de rotación calculada como el producto de las 3 matrices
    R = linalg.multi_dot([R_z_phi,R_y_tita,R_z_psi])
    #R = np.round(linalg.multi_dot([R_z_phi,R_y_tita,R_z_psi]),cant_de_decimales)
    
    # El índice de configuración
    #sg_tita = sign(tita)
    
    
    sg_tita = 1 if sign(tita) >= 0 else -1
    return R,sg_tita 




# In[]
"""
RMat2Eul es la función que resuelve el problema inverso de los ángulos de Euler
Los parámetros de entrada son la matriz de rotación R, el indice de configuración, y 
el ángulo phi actual
Ejemplo de como utilizar la función para una matriz de rotación R_rot, indice de
configuración, ind_conf y un ángulo actual ang_act
    ang_euler = RMat2Eul(R = R_rot, sg_tita = ind_conf, phi_act = ang_act )
"""
def RMat2Eul(R=eye(3),g = 1,phi_act = 0):
    
    # Defino todos los elementos de la matriz de rotación
    nx = R[0][0]; ny = R[1][0]; nz = R[2][0]
    sx = R[0][1]; sy = R[1][1]; sz = R[2][1]
    ax = R[0][2]; ay = R[1][2]; az = R[2][2]
    
    phi_act = phi_act*pi/180
    
    #phi_act = get_in_range(phi_act)
    # Si los valores ax y ay no son cero calculo las 2 posibles soluciones
    
    #if(ax != 0 or ay != 0):
    #if(round(ax,18) != 0 or round(ay,18) != 0):
    if( absolute(ax) > tol or absolute(ay) > tol): # Se va a usar esta línea en ves de != 0 para evitar futuros errores
        phi = arctan2(g*ay,g*ax)
    else:
        phi = phi_act
    
    #tita = arctan2(g*ax*cos(phi) + g*ay*sin(phi),az)
    #psi = arctan2( g*(-sin(phi)*nx + cos(phi)*ny) , g*( -sin(phi)*sx + cos(phi)*sy))
    
    tita = arctan2(ax*cos(phi) + ay*sin(phi),az)
    psi = arctan2( -sin(phi)*nx + cos(phi)*ny , -sin(phi)*sx + cos(phi)*sy)
    
    
   
    ang = np.array([phi,tita,psi])
    return ang*180/pi
        