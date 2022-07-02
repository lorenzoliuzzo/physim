
# author:          Lorenzo Liuzzo
# email:           lorenzoliuzzo@outlook.com
# description:     Different types of oscillators.
# last updated:    02/07/2022


from ode import ODE
import numpy as np


class HarmonicOscillator(ODE) : 

    def __init__(self, omega0) :
        self.omega0 = omega0

    def eval(self, pos, t = None) : 
        self.df[:, 0] = pos[:, 1]
        self.df[:, 1] = - np.dot(self.omega0 ** 2, pos[:, 0])
        return self.df

    def set_omega0(self, omega0) : 
        self.omega0 = omega0

    def get_omega0(self) : 
        return self.omega0 

    def print(self) : 
        print("omega0 = ", self.omega0)



class ForcedOscillator(ODE) : 

    def __init__(self, omega0, omega1) :
        self.omega0 = omega0
        self.omega1 = omega1

    def eval(self, pos, t) :
        self.df[:, 0] = pos[:, 1]
        self.df[:, 1] = - np.dot(self.omega0 ** 2, pos[:, 0]) + np.sin(self.omega1 * t)
        return self.df

    def set_omega0(self, omega0) : 
        self.omega0 = omega0

    def set_omega1(self, omega1) : 
        self.omega1 = omega1

    def get_omega0(self) : 
        return self.omega0 
        
    def get_omega1(self) : 
        return self.omega1 
        
    def print(self) : 
        print("omega0 = ", self.omega0)
        print("omega1 = ", self.omega1)



class DampedOscillator(ODE) : 
    def __init__(self, omega0, alpha) :
        self.omega0 = omega0
        self.alpha = alpha

    def eval(self, pos, t = None) : 
        self.df[:, 0] = pos[:, 1]
        self.df[:, 1] = - np.dot(self.omega0 ** 2, pos[:, 0]) - np.dot(self.alpha, pos[:, 1])
        return self.df

    def set_omega0(self, omega0) : 
        self.omega0 = omega0

    def set_alpha(self, alpha) : 
        self.alpha = alpha

    def get_omega0(self) : 
        return self.omega0 
        
    def get_alpha(self) : 
        return self.alpha 

    def print(self) : 
        print("omega0 = ", self.omega0)
        print("alpha = ", self.alpha)
        self.print_position()



class ForcedDampedOscillator(ODE) : 
    def __init__(self, omega0, omega1, alpha) :
        self.omega0 = omega0
        self.omega1 = omega1
        self.alpha = alpha

    def eval(self, pos, t) : 
        self.df[:, 0] = pos[:, 1]
        self.df[:, 1] = - np.dot(self.omega0 ** 2, pos[:, 0]) + np.sin(self.omega1 * t) - np.dot(self.alpha, pos[:, 1])
        return self.df

    def set_omega0(self, omega0) : 
        self.omega0 = omega0

    def set_omega1(self, omega1) : 
        self.omega1 = omega1

    def set_alpha(self, alpha) : 
        self.alpha = alpha

    def get_omega0(self) : 
        return self.omega0 
        
    def get_omega1(self) : 
        return self.omega1 
        
    def get_alpha(self) : 
        return self.alpha 

    def print(self) : 
        print("omega0 = ", self.omega0)

