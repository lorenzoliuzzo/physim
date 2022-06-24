
# author:          Lorenzo Liuzzo
# email:           lorenzoliuzzo@outlook.com
# description:     Different types of oscillators.
# last updated:    24/06/2022

from ode import ODE
import numpy as np


class HarmonicOscillator(ODE) : 

    def __init__(self, omega0) :
        self.omega0 = omega0

    def eval(self, pos = None, t = None) : 
        if pos is None : pos = self.get_position()
        self.df[:, 0] = pos[:, 1]
        self.df[:, 1] = - np.dot(pos[:, 0], self.omega0 ** 2)
        return self.df

    def print(self) : 
        print("omega0 = ", self.omega0)
        self.print_position()



# class ForcedOscillator(ODE) : 

#     def __init__(self, omega0, omega1) :
#         self.omega0 = omega0
#         self.omega1 = omega1

#     def eval(self, pos = None, t = None) :
#         if pos is None : pos = self.get_position()
#         self.df[:, 0] = pos[:, 1]
#         self.df[:, 1] = - np.dot(pos[:, 0], self.omega0 ** 2) - np.sin(self.omega1 * t)
#         return self.df

#     def print(self) : 
#         print("omega0 = ", self.omega0)
#         print("omega1 = ", self.omega1)
#         self.print_position()



class DampedOscillator(ODE) : 
    def __init__(self, omega0, alpha) :
        self.omega0 = omega0
        self.alpha = alpha


    def eval(self, pos = None, t = None) : 
        if pos is None : pos = self.get_position()
        self.df[:, 0] = pos[:, 1]
        self.df[:, 1] = - np.dot(self.omega0 ** 2, pos[:, 0]) + np.dot(self.alpha, pos[:, 1])
        return self.df

    def print(self) : 
        print("omega0 = ", self.omega0)
        print("alpha = ", self.alpha)
        self.print_position()



class ForcedDampedOscillator(ODE) : 
    def __init__(self, omega0, omega1, alpha) :
        self.omega0 = omega0
        self.omega1 = omega1
        self.alpha = alpha

    def eval(self, t, pos = None) : 
        if pos is None : pos = self.get_position()
        self.df[:, 0] = pos[:, 1]
        self.df[:, 1] = - np.dot(self.omega0 ** 2, pos[:, 0]) + np.sin(self.omega1 * t) - np.dot(self.alpha, pos[:, 1])
        return self.df

    def print(self) : 
        print("omega0 = ", self.omega0)
        self.print_position()

