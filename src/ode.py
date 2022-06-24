
# author:          Lorenzo Liuzzo
# email:           lorenzoliuzzo@outlook.com
# description:     Ordinaries Differential Equations.
# last updated:    24/06/2022

import numpy as np
from abc import ABC, abstractmethod
from position import Position 


class ODE(Position, ABC) :

    # differential
    df = np.empty(shape = (3, 2), dtype = float) 

    # __init__ with a position [coord, vel]
    @abstractmethod
    def __init__(self, coord, vel) :
        super().__init__()
        pass

    # evaluting the derivative of the position
    @abstractmethod
    def eval(self, t = None, pos = None) : 
        pass

    # Euler's method
    def euler(self, h = 0.001, t = None, pos = None) :
        pos = pos or self.get_position()
        appo = pos + h * self.eval(t, pos)
        self.set_position(appo[:, 0], appo[:, 1])
    
    # Runge Kutta's 4th order method
    def rk4(self, h = 0.001, t = None, pos = None) :
        t = t or 0
        pos = pos or self.get_position()
        k0 = np.empty(shape = (3, 2), dtype = float) 
        k1 = np.empty(shape = (3, 2), dtype = float) 
        k2 = np.empty(shape = (3, 2), dtype = float) 
        k3 = np.empty(shape = (3, 2), dtype = float) 
        
        k0 = self.eval(t, pos)
        k1 = self.eval(t + h / 2, pos + k0 * h / 2)
        k2 = self.eval(t + h / 2, pos + k1 * h / 2)
        k3 = self.eval(t + h, pos + k2 * h)   
        pos += h * (k0 + 2 * k1 + 2 * k2 + k3) / 6
        self.set_position(pos[:, 0], pos[:, 1])
            

