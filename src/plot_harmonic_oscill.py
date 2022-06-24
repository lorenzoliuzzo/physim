
# author:          Lorenzo Liuzzo
# email:           lorenzoliuzzo@outlook.com
# description:     Plot of an harmonic oscillator using ODE's numerical methods.
# last updated:    24/06/2022

import matplotlib.pyplot as plt
from oscillators import HarmonicOscillator


t = 0
tmax = 100

time = []
coord = []

print("This is a simulation of an harmonic oscillator.")

# generating an image
fig = plt.figure()
plt.xlabel('time [s]')
plt.ylabel('position [m]')
plt.title('Harmonic Oscillator')

# set oscillator
oscill = HarmonicOscillator(1)

# simulation with different increment 
for h in {0.0001} :

    # set position
    oscill.set_position([0, 0, 0], [1, 0, 0])

    while (t < tmax) :
        # updating data for plotting
        time.append(t)
        coord.append(oscill.get_coord_x())

        # ode solver
        oscill.rk4(h)
        t += h

    # plot position
    plt.plot(time, coord)
    
    t = 0
    time = []
    coord = []
    
print("Simulation ended!") 
print("Image generated!")
fig.savefig('HarmonicOscillator.jpg', bbox_inches='tight')
plt.show()
