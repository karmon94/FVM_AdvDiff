# -*- coding: utf-8 -*-
"""
Created on Mon Jun  4 13:47:58 2018

@author: karmon
"""

from mpi4py import MPI
import numpy as np 

LEFT, RIGHT = 0, 1
neighbour_processes = [0, 0]


class CartComm():
    
    """ Clase designada para construir una topologia cartesiana para ser usada
    en dominios simples rectangulares. Los subdominios de que genera esta topologia
    igualmente son rectangulares. La forma de la topologia es como sigue:
    
     ^ ^ eje - y (grid_rows)
     | |
     | +---------+---------+---------+---------+   
       |         |         |         |         |   
     J |    0    |    3    |    6    |    9    |   
       |  (0,0)  |  (1,0)  |  (2,0)  |  (3,0)  |   
     | |         |         |         |         |   
     | +---------+---------+---------+---------+---> eje - x  (grid_cols)
     V   
     """
    
    # Variables necesarias para programar con MPI #
    comm = MPI.COMM_WORLD
    rank = comm.rank
    size = comm.size

    # Variables designadas para la construccion de la topologia virtual #
    ##### Al ser 1D el valor de la fila tendra siempre valor de 1
    ##### Mientras que el numero de columnas dependera de los procesos ingresados por el usuario
    grid_rows = int(1)
    grid_cols = comm.size

    if rank == 0:
        print('Building a %d x %d grid topology' % (grid_rows, grid_cols))

    ##### Instruccion para crear la topologia virtual
    ##### (grid_rows, grid_cols) === numero de subdominios en que dividira el dominio original
    ##### periodo = activar o desactivar los ciclos en eje x o y
    cartesian_communicator = comm.Create_cart(
        (grid_rows, grid_cols), periods=(False, False), reorder=True)
    
    ##### Obtener las coordenadas del subdominio actual   (row, column)
    my_mpi_row, my_mpi_col = cartesian_communicator.Get_coords(
        cartesian_communicator.rank)
        
    ##### Comunicacion entre los subdominios adyacentes al subdominio actual
    ##### Si el valor obtenido es -1 no se tiene ningun subdominio
    neighbour_processes[LEFT], neighbour_processes[RIGHT] =\
        cartesian_communicator.Shift(1, 1)
        
    


if __name__ == '__main__':
    
    topology = CartComm()
    
    print('Process = %s \t \
          row = %s \t \
          column = %s \n \
          neighboor_processes[LEFT] = %s \n \
          neighboor_processes[RIGHT] = %s \n \
          neighboor = %s ' % (topology.rank, topology.my_mpi_row, topology.my_mpi_col,
          neighbour_processes[LEFT],
          neighbour_processes[RIGHT],
          neighbour_processes))
    