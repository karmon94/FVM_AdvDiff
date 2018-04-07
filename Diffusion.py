"""
Created on Tue Mar  6 15:11:05 2018

@author: Dr. Luis Miguel de la Cruz
adaptado por: Rogelio Carrillo

"""
# Importar Script donde se encuentra la parte
# de los Coeficienetes
from Coefficients import Coefficients

# Declaración de la Clase Difusion

class Diffusion1D(Coefficients):
    
    """
    Esta clase realiza el calculo en los arreglos para la parte Difusiva con una Gamma dada
    con los arreglos principales de los coeficientes del metodo de Volumen Finito.
    """ 
    ##### Declaaración de las Varibles de Numero y Tamaño de Volumenes
    # en conjunto con la variable Gamma para la parte Difusiva
    ### Inicialización####
    def __init__(self, nvx = None, Gamma = None, dx = None):
        """
        Constructor de la clase.
        
        @param nvx: número de volúmenes
        @param Gamma: coeficiente Gamma de la función a resolver
        @param dx: valor de los intervalos en x
        """
        super().__init__(nvx, dx)
        self.__nvx = nvx
        self.__Gamma = Gamma
        self.__dx = dx

    def __del__(self):
        """
        Destructor de la clase
        """
        del(self.__Gamma)
        del(self.__dx)
    
    # Calculo de los Coeficientes con la Gamma####
    def calcCoef(self):
        """
        Método que calcula la parte difusiva de los coeficientes que aparecerán en la matriz que se resolverá.
        """
        aE = self.aE()
        aW = self.aW()
        aP = self.aP()
                
        aE += self.__Gamma / self.__dx
        aW += self.__Gamma / self.__dx
        aP += aE + aW
        
        return aE, aW, aP



# Parte Principal del Script donde se invica la clase e Impresión de Pruebas
if __name__ == '__main__':
    
    df1 = Diffusion1D(5, 5, 1)
    df1.alloc(5)
    df1.calcCoef()
    df1.setSu(100)

    print('-' * 20)  
    print(df1.aP(), df1.aE(), df1.aW(), df1.Su(), sep = '\n')
    print('-' * 20)  

    df1.bcDirichlet('LEFT_WALL', 2)
    df1.bcDirichlet('RIGHT_WALL', 1)
    print(df1.aP(), df1.aE(), df1.aW(), df1.Su(), sep = '\n')
    print('-' * 20)  
