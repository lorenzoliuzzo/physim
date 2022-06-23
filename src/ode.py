
# author:          Lorenzo Liuzzo
# email:           lorenzoliuzzo@outlook.com
# description:     
# last updated:    23/06/2022

import numpy as np
from abc import ABC, abstractmethod
from position import Position 


class ODE(Position, ABC) :

    # differential
    df = np.empty(shape = (3, 2)) 

    # __init__ with a position [coord, vel]
    @abstractmethod
    def __init__(self, coord, vel) :
        super().__init__()
        pass

    # evaluting the derivative of the position
    @abstractmethod
    def eval(self, pos = None, t = None) : 
        pass

    # Euler's method
    def simpleEuler(self, h = 0.001, pos = None, t = None) :
        if (pos == None) : 
            appo = self.get_position() + h * self.eval(t)
            self.set_position(appo[:, 0], appo[:, 1])
        else : 
            appo = pos + h * self.eval(pos, t)
            self.set_position(appo[:, 0], appo[:, 1])
    
    # Runge Kutta's 4th order method
    def rk4(self, h = 0.001, pos = None, t = None) :
        if (pos == None) : 
            pos = self.get_position()

        k0 = self.eval(pos, t)
        k1 = self.eval(pos + k0 * h / 2,  t + h / 2)
        k2 = self.eval(pos + k1 * h / 2,  t + h / 2)
        k3 = self.eval(pos + k2 * h, t + h) 
        pos += h * (k0 + 2 * k1 + 2 * k2 + k3) / 6
        self.set_position(pos[:, 0], pos[:, 1])
            

