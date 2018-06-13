# -*- coding: utf-8 -*-
"""
Created on Wed Jun  6 16:25:01 2018

@author: karmon
"""

import FiniteVolumeMethod as fvm
#Librerias de Python
import numpy as np
import matplotlib.pyplot as plt
import CartComm as cart
from mpi4py import MPI
import time

# Variables necesarias para programar con MPI #
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


# Datos del Problema (Variables)
longitud = 0.5 # metros
TA = 100 # °C
TAa = 0
TB = 500 # °C 
TBa = 0
k  = 1000 # W/m.K
N1  = 50 # Número de nodos
N = N1
#
# Creamos la malla y obtenemos datos importantes
#
cart.CartComm()
    
dim = cart.CartComm.grid_cols
derecha = cart.RIGHT
izquierda = cart.LEFT


malla = fvm.Mesh(nodes = N, length = longitud)
nx    = malla.nodes()     # Número de nodos
nvx   = malla.volumes()   # Número de volúmenes
delta = malla.delta()     # Tamaño de los volúmenes
#
# Se construye un vector de coordenadas del dominio
#
x = malla.createMesh()
# 
###############################################################################
############# Si solo se selecciona un solo procesos ########################
###############################################################################
t1 = time.time()
if dim == 1:
    malla = fvm.Mesh(nodes = N, length = longitud)
    nx    = malla.nodes()     # Número de nodos
    nvx   = malla.volumes()   # Número de volúmenes
    delta = malla.delta()     # Tamaño de los volúmenes

    #
    #  Creamos los coeficientes del Metodo de Volumen Finito
    #
    df1 = fvm.Diffusion1D(nvx, Gamma = k, dx = delta) #Invocacion del Metodo de Volumen Finito
    df1.alloc(nvx) # Se aloja memoria para los coeficientes
    df1.calcCoef() # Se calculan los coeficientes
    
    

    data = T = np.zeros(nvx) # El arreglo contiene ceros
    print ("Solucionaremos la T = ", data, "en el rank =", rank)
    T[0]  = TA        # Condición de frontera izquierda
    T[-1] = TB       # Condición de frontera derecha
    df1.bcDirichlet('LEFT_WALL', T[0])   # Se actualizan los coeficientes
    df1.bcDirichlet('RIGHT_WALL', T[-1]) # de acuerdo a las cond. de frontera
    print ("ahora t vale = ", T)
    
    Su = df1.Su()  # Vector del lado derecho
    A = fvm.Matrix(malla.volumes())  # Matriz del sistema
    A.build(df1) # Construcción de la matriz en la memoria
    # Se resuelve el sistema usando un algoritmo del módulo linalg
    T[1:-1] = np.linalg.solve(A.mat(),Su[1:-1])
    # Impresión del vector Solución de la Matriz construida
    print('Solución para T= {}'.format(T))
    print('.'+'-'*70+'.')
    t2 = time.time()
    print("Tiempo realizado: ")
    print('MPI: ' + str((t2-t1)) + "\n")
    #
    # Se construye un vector de coordenadas del dominio
    #
    x = malla.createMesh()
    #
    # Calculamos la solución analítica
    #
    Ta = 800 * x + 100
    #
    #  Se grafica la solución
    #
    x *= 100 # Transformación a [cm]
    plt.plot(x,Ta, '-', label = 'Sol. analítica') # Sol. analítica
    plt.plot(x,T,'o', label = 'Sol. FVM')
    plt.title('Solución de $k (\partial^2 T/\partial x^2) = 0$ con FVM usando MPI')
    plt.xlabel('$x$ [cm]')
    plt.ylabel('Temperatura [$^o$C]')
    plt.grid()
    plt.legend()
    plt.savefig('Tarea_ejemplo4.1.pdf')
    plt.show()

###############################################################################
################ Si se seleccionan dos procesos ###########################
###############################################################################
elif dim == 2:
    
    N  = int((N1/2)+1)
    
    t1 = time.time()
    
    malla = fvm.Mesh(nodes = N, length = longitud)
    nx    = malla.nodes()     # Número de nodos
    nvx   = malla.volumes()   # Número de volúmenes
    delta = malla.delta()     # Tamaño de los volúmenes
    
    
    #
    #  Creamos los coeficientes del Metodo de Volumen Finito
    #
    df1 = fvm.Diffusion1D(nvx, Gamma = k, dx = delta) #Invocacion del Metodo de Volumen Finito
    df1.alloc(nvx) # Se aloja memoria para los coeficientes
    df1.calcCoef() # Se calculan los coeficientes
    
    if rank == 0:
        data = T = np.zeros(nvx) # El arreglo contiene ceros
        print ("Solucionaremos la T = ", data, "en el rank =", rank)
        T[0]  = TA        # Condición de frontera izquierda
        T[-1] = TAa        # Condición de frontera derecha
        df1.bcDirichlet('LEFT_WALL', T[0])   # Se actualizan los coeficientes
        df1.bcDirichlet('RIGHT_WALL', T[-1]) # de acuerdo a las cond. de frontera
        print ("ahora t vale = ", T)
        
        Su = df1.Su()  # Vector del lado derecho
        
        A = fvm.Matrix((malla.volumes()))  # Matriz del sistema
        A.build(df1) # Construcción de la matriz en la memoria
        # Se resuelve el sistema usando un algoritmo del módulo linalg
        T[1:-1] = np.linalg.solve(A.mat(),Su[1:-1])
        # Impresión del vector Solución de la Matriz construida
        print('Solución para T= {}'.format(T))
        print('.'+'-'*70+'.')
        
        #######################################################################
        ####### Entre mas nodos, por ende mayor numero de iteraciones #########
        #######################################################################
        for i in range (200):
            #### enviando frontera izquierda fantasma a derecha fantasma ####  
            #########    rank = 0    real  ||--+---+---+---** f   --->>      rank = 1    f **--+---+---+---|| real
            dataizq = T[-2]
            comm.send(dataizq, dest = derecha, tag=11)
             
            #### recibiendo frontera derecha fantasma y asignando  a izquierda fantasma ###
            #########   rank = 1    f **--+---+---+---|| real   --->>    rank = 0    real  ||--+---+---+---** f     
            datader = comm.recv(source = derecha, tag =12)
                    
            df1.cleanCoefficients()
            df1.calcCoef()
             
            T[-1]  = datader# Condición de frontera izquierda
            T[0] = TA        # Condición de frontera derecha
            df1.bcDirichlet('LEFT_WALL', T[0])   # Se actualizan los coeficientes
            df1.bcDirichlet('RIGHT_WALL', T[-1]) # de acuerdo a las cond. de frontera
            Su = df1.Su()  # Vector del lado derecho
            A = fvm.Matrix(malla.volumes())  # Matriz del sistema
            A.build(df1) # Construcción de la matriz en la memoria
            T[1:-1] = np.linalg.solve(A.mat(),Su[1:-1])
            # Impresión del vector Solución de la Matriz construida
            #print('Solución para T = {}'.format(T))
            #print('.'+'-'*70+'.')
    
        ################ recibiendo el arreglo final del rank = 1 #############   
        dataTotal = comm.recv(source = 1, tag =13)
        
        #######################################################################
        ######## Se juntan los arreglos sin las fronteras fantasmas  ##########
        #######################################################################
        TF = T[:-2]
        TTF = dataTotal[1:]
        Solucion = []
        Solucion.extend(TF)
        Solucion.extend(TTF)
        print ("Solucion = ", Solucion)
        
        #######################################################################
        ##### Obtencion del tiempo de ejecucion del obtencion de solucion #####
        #######################################################################
        t2 = time.time()
        print("Tiempo realizado: ")
        print('MPI: ' + str((t2-t1)) + "\n")
        
        
        # Calculamos la solución analítica
        #
        Ta = 800 * x + 100
        #
        #  Se grafica la solución
        #
        x *= 100 # Transformación a [cm]
        plt.plot(x,Ta, '-', label = 'Sol. analítica') # Sol. analítica
        plt.plot(x,Solucion,'o', label = 'Sol. FVM')
        plt.title('Solución de $k (\partial^2 T/\partial x^2) = 0$ con FVM usando MPI')
        plt.xlabel('$x$ [cm]')
        plt.ylabel('Temperatura [$^o$C]')
        plt.grid()
        plt.legend()
        plt.savefig('Tarea_ejemplo4.1.pdf')
        plt.show()
        
###############################################################################
##########################   RANK 2   #########################################        
###############################################################################        
        
    else:
        
        data = TT = np.zeros(nvx) # El arreglo contiene ceros
        print ("Solucionaremos la T = ", data, "en el rank =", rank)
        TT[0]  = TBa        # Condición de frontera izquierda
        TT[-1] = TB        # Condición de frontera derecha
        df1.bcDirichlet('LEFT_WALL', TT[0])   # Se actualizan los coeficientes
        df1.bcDirichlet('RIGHT_WALL', TT[-1]) # de acuerdo a las cond. de frontera
        print ("ahora TT vale = ", TT)
        
        Su = df1.Su()  # Vector del lado derecho
        A = fvm.Matrix(malla.volumes())  # Matriz del sistema
        A.build(df1) # Construcción de la matriz en la memoria
        # Se resuelve el sistema usando un algoritmo del módulo linalg
        TT[1:-1] = np.linalg.solve(A.mat(),Su[1:-1])
        # Impresión del vector Solución de la Matriz construida
        print('Solución para TT= {}'.format(TT))
        print('.'+'-'*70+'.')
        
        #######################################################################
        ####### Entre mas nodos, por ende mayor numero de iteraciones #########
        #######################################################################
        for i in range (200):
            ######## enviando frontera fantasma lado derecha a frontera lado izquierdo #########     
            #########    rank = 1    f **--+---+---+---|| real   --->     rank = 0    real  ||--+---+---+---** f       
            datader = TT[1]
            comm.send(datader, dest=0, tag=12)
             
            ##### recibiendo frontera izquierda fantasma y asignando frontera derecha fantasma #####
            #########    rank = 0    real  ||--+---+---+---** f   --->>      rank = 1    f **--+---+---+---|| real
            dataizq = comm.recv (source = 0, tag =11)
            
            df1.cleanCoefficients()
            df1.calcCoef()
             
            TT[0]  = dataizq        # Condición de frontera izquierda
            TT[-1] = TB        # Condición de frontera derecha
            df1.bcDirichlet('LEFT_WALL', TT[0])   # Se actualizan los coeficientes
            df1.bcDirichlet('RIGHT_WALL', TT[-1]) # de acuerdo a las cond. de frontera
            Su = df1.Su()  # Vector del lado derecho
            A = fvm.Matrix(malla.volumes())  # Matriz del sistema
            A.build(df1) # Construcción de la matriz en la memoria
            TT[1:-1] = np.linalg.solve(A.mat(),Su[1:-1])
            # Impresión del vector Solución de la Matriz construida
            #print('Solución para TT = {}'.format(TT))
            #print('.'+'-'*70+'.')
        
        #######################################################################
        ############# Envio final de la solucion al rank principal ############
        #######################################################################
        dataTotal = TT
        comm.send(dataTotal, dest=0, tag=13)
###############################################################################

