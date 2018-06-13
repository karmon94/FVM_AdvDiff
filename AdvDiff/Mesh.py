"""
Created on Tue Mar  6 15:11:05 2018

@author: Dr. Luis Miguel de la Cruz
adaptado por: Rogelio Carrillo

"""
### Librerias de Python ###
import numpy as np

class Mesh():
    
    """
    Esta clase realiza la construcción y/o ajuste de la malla para 
    calcular la solución del sistema, apartir de los 
    de los datos ingresados en el sistema,como numero de nodos
     y la longitud, volumenes, etc.
    """
    #### Inicialización de cada uno de las variables y funciones ##########
    # Se realiza de la siguiente manera con el fin de realizar
    # Con Programación Orientada a Objetos
    def __init__(self, nodes = None, 
                     volumes = None,
                     length = None):
        """
        Constructor de la clase.
        
        @param nodes: número de nodos
        @param volumes: número de volumenes
        @param lenght: longitud de la malla
        """
        
        self.__nodes = nodes
        self.__volumes = volumes
        self.__length = length     
        self.__delta = 1
        self.adjustNodesVolumes(nodes, volumes)
        self.calcDelta()
    
    def __del__(self):
        """
        Destructor de la clase
        """
        del(self.__nodes)
        del(self.__volumes)
        del(self.__length)
        del(self.__delta)
        
    def adjustNodesVolumes(self,nodes,volumes):
        """
        Metodo que ajusta el numero de nodos o volumenes, dependiendo de lo que ingrese el usuario
        """
        if nodes:
            self.__volumes = self.__nodes + 1
        if volumes:
            self.__nodes = self.__volumes - 1        
        
    def nodes(self):
        """ Metodo que regresa el numero de nodos
        """
        return self.__nodes
    
    def setNodes(self, nodes):
        """Metodo que asigna el numero de nodos
        """
        self.__nodes = nodes
        self.adjustNodesVolumes(nodes = nodes, volumes = None)
        
    def volumes(self):
        """ Metodo que regresa el numero de volumenes
        """
        return self.__volumes

    def setVolumes(self, volumes):
        """ Metodo que asigna el numero de nodos
        """
        self.__volumes = volumes
        self.adjustNodesVolumes(nodes = None, volumes = volumes)
        
    def length(self):
        """ Metodo que regresa el valor de la longitud
        """
        return self.__length
        
    def calcDelta(self):
        """ Metodo que calcula el valor de los intervalos dependiendo de la longitud y los nodos 
        """
        if self.__length:
            self.__delta = self.__length / (self.__nodes - 1)
        
    def delta(self):
        """ Metodo que regresa el valor de los intervalos espaciales
        """
        return self.__delta
    #################################################
    
    ### Creación de la Malla####
    def createMesh(self):
        """ Metodo que construye la malla
        """
        first_volume = self.__delta / 2
        final_volume = self.__length - first_volume
        self.__x = np.zeros(self.__volumes)
        self.__x[1:-1] = np.linspace(first_volume,final_volume,self.__volumes-2)
        self.__x[-1] = self.__length
        return self.__x

# Parte Principal del Script donde se invica la clase e Impresión de Pruebas de donde estan los nodos y volumenes  
if __name__ == '__main__':

    m1 = Mesh()
    print(m1.nodes(), m1.volumes())
    print('_' * 20)   
    
    m1 = Mesh(nodes = 5)
    print(m1.nodes(), m1.volumes())
    print('_' * 20)
   
    m1 = Mesh(volumes = 5)
    print(m1.nodes(), m1.volumes())
    print('_' * 20)
    
    m1 = Mesh(5,5)
    print(m1.nodes(), m1.volumes())
    print('_' * 20)
    
    m1.setNodes(8)
    print(m1.nodes(), m1.volumes())
    print('_' * 20)

    m1.setVolumes(8)
    print(m1.nodes(), m1.volumes())
    print('_' * 20)
    
    m1 = Mesh(nodes =  5, length = 33)
    print(m1.nodes(), m1.volumes(), m1.length())
    print('_' * 20)
    
    m1 = Mesh(volumes =  5, length = 33)
    print(m1.nodes(), m1.volumes(), m1.length())
    print('_' * 20)
    
    m1 = Mesh(nodes = 5, length = 1)
    print(m1.nodes(), m1.volumes(), m1.length(), m1.delta())
    print('_' * 20)    
    
    m1 = Mesh(volumes = 10, length = 1)
    print(m1.nodes(), m1.volumes(), m1.length(), m1.delta())
    print('_' * 20) 

    m1 = Mesh(volumes = 6, length = 1)
    print(m1.nodes(), m1.volumes(), m1.length(), m1.delta())
    m1.createMesh()
    print('_' * 20) 
    