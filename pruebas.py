#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 16:04:35 2019

@author: franco
"""

# In[]


# Se importa la biblioteca numpy para aplicar operaciones similares a las que se usan en matlab
# Adenás esta biblioteca tiene funciones para cargar datos desde un archivo csv
import numpy as np # importa numpy y se puede acceder a la biblioteca con el apodo np
from tp2 import * # Importa todas las funciones del archivo 'tp1.py'
# In[]

# Cargo el archivo de prueba con la función loadtxt
file_name = "variables_articulares_prueba.txt"
datos_prueba = np.loadtxt(file_name, delimiter=',', skiprows=1)

# Verificación: Primero se calcula la matrix de rotación y luego se aplica la operación inversa
for ind, var_art in enumerate(datos_prueba):
    t1_act = var_art[0] # Defino phi_actu como el valor de phi que se están probando
    t4_act = var_art[3]
    A, g, t_act= pos_prob_dir(*var_art) # Calculo la matriz de rotación y el índice de configuración
    
    """
    print("DEBUG")
    print("Línea " + str(ind+2))
    print(var_art)
    print(A)
    print(g)
    print(t_act)
    print("\n\n")
    """
    # Verifico si al resolver el problema inverso, se obtiene los mismos ángulos
    
    # Calculo los ángulos de euler a partir de la matriz de rotación obtenida
    var_art_res = pos_prob_inv(A,*g,*t_act)
    #print(var_art_res)
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
    
    if ( not np.all(np.abs(var_art_res - var_art) < tol)):
        print("Error para los ángulos que se encuentran en la línea: " + str(ind + 2))