#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 08:39:11 2019

@author: franco
"""

# In[]
# Se importa la biblioteca numpy para aplicar operaciones similares a las que se usan en matlab
from numpy import *

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
    
    # Se escriben las matrices de rotación y con la función array se convierten en elementos
    # de numpy
    R_z_phi = [[cos(phi), -sin(phi), 0],[sin(phi), cos(phi), 0],[0,0,1]]
    R_z_phi = array(R_z_phi)
    #print(R_z_phi)
    
    R_y_tita = [[cos(tita), 0, sin(tita)], [0,1,0] , [-sin(tita), 0, cos(tita)]]
    R_y_tita = array(R_y_tita)
    #print(R_y_tita)
    
    R_z_psi = [[cos(psi), -sin(psi), 0],[sin(psi), cos(psi), 0],[0,0,1]]
    R_z_psi = array(R_z_psi)
    #print(R_z_psi)
    
    # Matriz de rotación calculada como el producto de las 3 matrices
    R = linalg.multi_dot([R_z_phi,R_y_tita,R_z_psi])
    # El índice de configuración
    sg_tita = sign(tita)
    
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
def RMat2Eul(R=eye(3),sg_tita = 1,phi_act = 0):
    
    # Defino todos los elementos de la matriz de rotación
    nx = R[0][0]; ny = R[1][0]; nz = R[2][0]
    sx = R[0][1]; sy = R[1][1]; sz = R[2][1]
    ax = R[0][2]; ay = R[1][2]; az = R[2][2]
    
    tol = 1e-12
    
    # Si los valores ax y ay no son cero calculo las 2 posibles soluciones
    #if(ax != 0 or ay != 0):
    if( absolute(ax) > tol or absolute(ay) > tol): # Se va a usar esta línea en ves de != 0 para evitar futuros errores
        # Calculo las 2 posibles soluciones de phi teniendo en cuenta que 
        # arctan2 devuelve valores entre -pi y pi y que los valores de phi también
        # están definidos en ese rango                
        phi1 = arctan2(ay,ax)
        # Se agrega la otra posible solución dependiendo si la primera es
        # negativa o positiva 
        phi2 = (phi1 + pi) if phi1<0 else (phi1 - pi) 
        phi = array([phi1, phi2])  # Se arma un arreglo con ambas soluciones
        
        # Calculo ambos valores posibles de tita
        tita1 = arctan2(ax*cos(phi1) + ay*sin(phi1),az)
        tita2 = -tita1
        tita = array([tita1,tita2])
        
        # Calculo ambos valores posibles de psi
        psi1 = arctan2(-sin(phi1)*nx + cos(phi1)*ny ,-sin(phi1)*sx + cos(phi1)*sy)
        psi2 = (psi1 + pi) if psi1<0 else (psi1 - pi)
        psi = array([psi1, psi2])
        
        # Si el indice de configuración es igual al signo de tita1
        # defino un vector con las soluciones 1, en caso contrario con  las soluciones 2
        ang = array([phi1,tita1,psi1]) if (sg_tita == sign(tita1) )  else array([phi2,tita2,psi2])
    else:
        # Si ax = ay = 0, entonces calculo el valor de psi a partir del valor de phi ingresados
        phi_act = phi_act*pi/180 
        psi = arctan2(-sin(phi_act)*nx + cos(phi_act)*ny ,-sin(phi_act)*sx + cos(phi_act)*sy)
        ang = array([phi_act,0,psi])
    
    return ang*180/pi 
        