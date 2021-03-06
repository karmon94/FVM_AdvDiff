#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  4 15:32:38 2018

@author: Rogelio Carrillo

Solución del Problema 5.1 del Libro "An Introduction to Computational Fluid Dynamics"
de H K Versteeg and W Malalasekera con Metodo Upwind 
----------------------------------


___________________________________________________________

 |<-------------- 1.0 m ------------->|
                                     
            u ---> 
 |------------------------------------|
                                 
\phi_0 = 1                       \phi_L = 0
___________________________________________________________



"""

#Importar Script donde se encuentra el Metodo de Volumen Finito
import FiniteVolumeMethod as fvm
#Librerias de Python
import numpy as np
import matplotlib.pyplot as plt

### Calculo Mediante la Solución Analitica
def analyticSol(x):
    return (np.exp(rho * u * x / Gamma) - 1) / (np.exp(rho * u * L / Gamma) - 1) * (phiB - phiA) + phiA

# Datos del Problema (Variables)
L = 1.0 # m
rho = 1.0 # kg/m^3
u = 0.1 # m/s
Gamma = 0.1 # kg / m.s
phiA = 1 #
phiB = 0 #
N = 6 # Número de nodos
#
# Creamos la malla y obtenemos datos importantes
#
malla = fvm.Mesh(nodes = N, length = L)
nx    = malla.nodes()     # Número de nodos
nvx   = malla.volumes()   # Número de volúmenes
delta = malla.delta()     # Tamaño de los volúmenes
#
# Imprimimos los datos del problema (nicely)
#
fvm.printData(Longitud = L,
              Densidad = rho,
              Velocidad = u,
              Coef_Diff = Gamma,
              Prop_0 = phiA,
              Prop_L = phiB,
              Nodos = nx, 
              Volúmenes = nvx,
              Delta = delta)
#
# Se aloja memoria para los coeficientes
#
coef = fvm.Coefficients()
coef.alloc(nvx)
#fvm.Coefficients.alloc(nvx)
#
#  Calculamos los coeficientes de FVM de la Difusión
#
dif = fvm.Diffusion1D(nvx, Gamma = Gamma, dx = delta)
dif.calcCoef()
#print('aW = {}'.format(dif.aW()), 
#      'aE = {}'.format(dif.aE()), 
#      'Su = {}'.format(dif.Su()), 
#      'aP = {}'.format(dif.aP()), sep='\n')
#print('.'+'-'*70+'.')
#
#  Calculamos los coeficientes de FVM de la Advección
#
adv = fvm.Advection1D(nvx, rho = rho, dx = delta)
adv.setU(u)
adv.calcCoef('Upwind') 
#print('aW = {}'.format(dif.aW()), 
#      'aE = {}'.format(dif.aE()), 
#      'Su = {}'.format(dif.Su()), 
#      'aP = {}'.format(dif.aP()), sep='\n')
#print('u = {}'.format(adv.u()))
#print('.'+'-'*70+'.')
#
# Se construye el arreglo donde se guardará la solución
#
phi = np.zeros(nvx) # El arreglo contiene ceros
phi[0]  = phiA       # Condición de frontera izquierda
phi[-1] = phiB       # Condición de frontera derecha
#
# Se aplican las condiciones de frontera
#
coef.bcDirichlet('LEFT_WALL', phi[0])   # Se actualizan los coeficientes
coef.bcDirichlet('RIGHT_WALL', phi[-1]) # de acuerdo a las cond. de frontera
#print('aW = {}'.format(dif.aW()), 
#      'aE = {}'.format(dif.aE()), 
#      'Su = {}'.format(dif.Su()), 
#      'aP = {}'.format(dif.aP()), sep='\n')
#print('u = {}'.format(adv.u()))
#print('.'+'-'*70+'.')
#
# Se construye el sistema lineal de ecuaciones a partir de los coef. de FVM
#
Su = coef.Su()  # Vector del lado derecho
A = fvm.Matrix(malla.volumes())  # Matriz del sistema
A.build(coef) # Construcción de la matriz en la memoria
print('A = ', A.mat(),
      'b = {}'.format(Su[1:-1]), sep='\n')
print('.'+'-'*70+'.')
# Se resuelve el sistema usando un algoritmo del módulo linalg
#
phi[1:-1] = np.linalg.solve(A.mat(),Su[1:-1])
print('Solución = {}'.format(phi))
print('.'+'-'*70+'.')
#
# Se construye un vector de coordenadas del dominio
#
x = malla.createMesh()
#
# Calculamos la solución exacta y el error
#
phi_a = analyticSol(x)
error = fvm.calcError(phi, phi_a)
datos = {'x(m)': x,
         'phi(x)': phi,
         'Analytic': phi_a,
         'Error': error}
fvm.printFrame(datos)
print('||Error|| = ', np.linalg.norm(error))
print('.'+ '-'*70 + '.')
#
#  Se grafica la solución
#
plt.plot(x,phi_a, '-', label = 'Sol. analítica') 
plt.plot(x,phi,'o', label = 'Sol. FVM')
plt.title('Solución de $\partial(p u \phi)/\partial x= \partial (\Gamma \partial\phi/\partial x)/\partial x$ con FVM')
plt.xlabel('$x$ [m]')
plt.ylabel('$\phi$ [...]')
plt.grid()
plt.legend()
plt.savefig('Tarea_ejemplo5.1(Upwind).pdf')
plt.show()
