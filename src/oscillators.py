from ode import ODE

class HarmonicOscillator(ODE) : 

    def __init__(self, omega0) :
        self.omega0 = omega0

    def eval(self, pos, t = None) : 
        self.df[:, 0] = self.get_velocity()
        self.df[:, 1] = np.dot(self.get_coordinates(), - (self.omega0 ** 2))
        return self.df

    def print(self) : 
        print("omega0 = ", self.omega0)
        self.print_position()


class ForcedOscillator(ODE) : 

    def __init__(self, omega0, omega1) :
        self.omega0 = omega0
        self.omega1 = omega1

    def eval(self, t) : 
        self.df[:, 0] = self.get_velocity()
        appo = - np.dot(self.get_coordinates(), (self.omega0 ** 2))
        self.df[:, 1] = appo + np.sin(self.omega1 * t)
        return self.df

    def print(self) : 
        print("omega0 = ", self.omega0)
        print("omega1 = ", self.omega1)
        self.print_position()
        

class DampedOscillator(ODE) : 
    def __init__(self, omega0, alpha) :
        self.omega0 = omega0
        self.alpha = alpha

    def eval(self, t = None) : 
        self.df[:, 0] = self.get_velocity()
        self.df[:, 1] = - (self.omega0 ** 2) * self.get_coordinates() - self.alpha * self.vel.get_velocity()
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

    def eval(self, t) : 
        self.df[:, 0] = self.get_velocity()
        self.df[:, 1] = - (self.omega0 ** 2) * self.get_coordinates() + np.sin(self.omega1 * t) - self.alpha * self.vel.get_velocity()
        return self.df

    def print(self) : 
        print("omega0 = ", self.omega0)
        self.print_position()

