
# author:          Lorenzo Liuzzo
# email:           lorenzoliuzzo@outlook.com
# description:     Ordinaries Differential Equations.
# last updated:    02/07/2022


import numpy as np
from abc import ABC, abstractmethod


class ODE(ABC) :

    # differential
    df = np.empty(shape = (3, 2), dtype = float) 

    # __init__ with a position [coord, vel]
    @abstractmethod
    def __init__(self) :
        pass

    # evaluting the derivative of the position
    @abstractmethod
    def eval(self, pos = None, t = None) : 
        pass

    # Euler's method
    def euler(self, pos, t = None, h = 0.001) :
        return pos.get_position() + h * self.eval(pos.get_position(), t)
    
    # Runge Kutta's 4th order method
    def rk4(self, pos, t = None, h = 0.001) :
        t = t or 0
        k0 = np.empty(shape = (3, 2), dtype = float)  
        k1 = np.empty(shape = (3, 2), dtype = float) 
        k2 = np.empty(shape = (3, 2), dtype = float) 
        k3 = np.empty(shape = (3, 2), dtype = float) 
        k0 = self.eval(pos.get_position(), t)
        k1 = self.eval(pos.get_position() + k0 * h / 2, t + h / 2)
        k2 = self.eval(pos.get_position() + k1 * h / 2, t + h / 2)
        k3 = self.eval(pos.get_position() + k2 * h, t + h)   
        return pos.get_position() + h * (k0 + 2 * k1 + 2 * k2 + k3) / 6

