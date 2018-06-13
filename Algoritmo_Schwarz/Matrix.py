"""
Created on Tue Mar  6 15:11:05 2018

@author: Dr. Luis Miguel de la Cruz
adaptado por: Rogelio Carrillo

"""
#Libreria de Python
import numpy as np

class Matrix():
    
    """
    Esta clase realiza la construcción de la matriz para 
    calcular la solución del sistema, apartir de los 
    coeficientes ya previamente calculados.
    """
    #### Inicialización de cada uno de las variables y funciones ##########
    # Se realiza de la siguiente manera con el fin de realizar
    # Con Programación Orientada a Objetos
    
    def __init__(self, nvx = None):
        """"Constructor de la clase
        """
        self.__N = nvx - 2 
        self.__A = np.eye(self.__N)

    def __del__(self):
        """Destructor de la clase
        """
        del(self.__N)
        del(self.__A)
        
    def mat(self):
        """ Metodo que regresa la matriz
        """
        return self.__A
    
    def build(self, coefficients = None):
        
        """ Metodo que construye la matriz con base a los coeficientes de un problema"""

        
# Esquema General de Como Estarioa Construida la Matriz en el caso de que
# nx = 5 y nvx = 6
# 0     1     2     3     4     5  <-- Volumes 
# o--|--x--|--x--|--x--|--x--|--o
#       0     1     2     3        <-- Unknowns    
#
#        0   1   2   3
# --+---------------------
# 0 | [[12. -4.  0.  0.]
# 1 | [ -4.  8. -4.  0.]
# 2 | [  0. -4.  8. -4.]
# 3 | [  0.  0. -4. 12.]]

        # Coeficientes##
        aP = coefficients.aP()
        aE = coefficients.aE()
        aW = coefficients.aW()
        A = self.__A
        
        #### Contrusción de la Matriz con Respecto a los Coeficientes####
        A[0][0] = aP[1]
        A[0][1] = -aE[1]
        for i in range(1,self.__N-1): # range(1,N-3)  <-- (1,2)
            A[i][i] = aP[i+1]
            A[i][i+1] = -aE[i+1]
            A[i][i-1] = -aW[i+1]
        A[-1][-1] = aP[-2]
        A[-1][-2] = -aW[-2]
        
        
        
        
# Parte Principal del Script donde se invica la clase e Impresión de Pruebas
if __name__ == '__main__':

    a = Matrix(6)
    print('-' * 20)  
    print(a.mat())
    print('-' * 20)  
    
    from Diffusion import Diffusion1D        
    df1 = Diffusion1D(6, 1, 0.25)
    df1.alloc(6)
    df1.calcCoef()
    df1.setSu(100)
    print(df1.aP(), df1.aE(), df1.aW(), df1.Su(), sep = '\n')
    print('-' * 20)  

    df1.bcDirichlet('LEFT_WALL', 2)
    df1.bcDirichlet('RIGHT_WALL', 1)
    print(df1.aP(), df1.aE(), df1.aW(), df1.Su(), sep = '\n')
    print('-' * 20)  
    
    a.build(df1)
    print(a.mat())
    print('-' * 20)  
