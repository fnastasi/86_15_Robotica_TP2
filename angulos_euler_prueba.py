#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 11:23:54 2019

@author: franco
"""

# In[]

# Se importa la biblioteca numpy para aplicar operaciones similares a las que se usan en matlab
# Adenás esta biblioteca tiene funciones para cargar datos desde un archivo csv
import numpy as np # importa numpy y se puede acceder a la biblioteca con el apodo np
from tp1 import * # Importa todas las funciones del archivo 'tp1.py'
# In[]

# Cargo el archivo de prueba con la función loadtxt
file_name = "angulos_prueba.txt"
datos_prueba = np.loadtxt(file_name, delimiter=',', skiprows=1)

# Verificación: Primero se calcula la matrix de rotación y luego se aplica la operación inversa
for ind, Euler_ang in enumerate(datos_prueba):
    phi_act = Euler_ang[0] # Defino phi_actu como el valor de phi que se están probando
    R, ind_conf = Eul2RMat(*Euler_ang) # Calculo la matriz de rotación y el índice de configuración
    
    """
    print("DEBUG")
    print("Línea " + str(ind+2))
    print(Euler_ang)
    print(RMat2Eul(R,ind_conf,phi_act))
    print("\n\n")
    """
    # Verifico si al resolver el problema inverso, se obtiene los mismos ángulos
    
    # Calculo los ángulos de euler a partir de la matriz de rotación obtenida
    Euler_ang_res = RMat2Eul(R,ind_conf,phi_act)
    
    tol = 1e-6 # Tolerancia
    
    
    
    # Comparo si al hacer la resta entre los valores de los ángulos originales del archivo
    # y lo que devuelve la función que resuelve el problema inverso
    # son menores en valor absoluto para un nivel de tolerancia.
    
    # Nota:
    # Euler_ang_res - Euler_ang devuelve un vector de 3 elementos de tipo float
    
    # np.abs(Euler_ang_res - Euler_ang) devuelve un vector de 3 elementos de tipo float 
    #en donde se le aplica la función módulo a cada elemento
    
    # np.abs(Euler_ang_res - Euler_ang) < tol devuelve un vector de 3 elementos de
    # de tipo booleano con True si se cumple la condición y False si no
    
    # p.all(np.abs(Euler_ang_res - Euler_ang) < tol) devuelve un valor booleano
    # Devuelve True si y solo si todos los elementos de  np.abs(Euler_ang_res - Euler_ang) < tol
    # son True
    
    if ( not np.all(np.abs(Euler_ang_res - Euler_ang) < tol)):
        print("Error para los ángulos que se encuentran en la línea: " + str(ind + 2))
        

# In[]

# Este script es para comparar las funciones que se escribieron en este tp 
# con el resultado de las funciones que tiene la biblioteca numpy
# para manejar rotaciones

### DOCUMENTACIÓN ####
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.transform.Rotation.html
        
# Se importa la biblioteca: Mediante la palabra clave 'as', podemos referirnos a la  bilioteca
# como Rot
from scipy.spatial.transform import Rotation as Rot

       
# Cargo el archivo de prueba con la función loadtxt
file_name = "angulos_prueba.txt"
datos_prueba = np.loadtxt(file_name, delimiter=',', skiprows=1) 
        


# Para cada conjunto de 3 ángulos de Euler en el archivo, comparo el resultado
# de la función que resuelve el problema directo con lo que devuelve la función
# de la biblioteca numpy
for ind, Euler_ang in enumerate(datos_prueba): # En Eulern_ang se guardan los 3 valores de los ángulos en cada interación
    
    # r es un objeto que contiene la información de la rotación, utilizando los ángulos de la variable
    # Euler_ang
    r = Rot.from_euler('ZYZ', Euler_ang, degrees=True)
    
    # Calculo la matriz de rotación utilizando la función que escribí
    R, ind_conf = Eul2RMat(*Euler_ang)
    
    tol = 1e-6 # Tolerancia
    
    # Comparo si al hacer la resta entre los valores de los elementos de la matriz de rotación
    # y lo que devuelve el objeto r, el resultado es menor a la tolerancia en valor absoluto
    
    
    if(not np.all(np.abs( R - r.as_dcm() ) < tol)):
        print("Error para los ángulos que se encuentran en la línea: " + str(ind + 2))
