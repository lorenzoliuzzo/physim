
# author:          Lorenzo Liuzzo
# email:           lorenzoliuzzo@outlook.com
# description:     Plot of a forced oscillator using ODE's numerical methods.
# last updated:    24/06/2022

import matplotlib.pyplot as plt
from oscillators import ForcedOscillator


t = 0
tmax = 100

time = []
coord = []

print("This is a simulation of a forced oscillator.")

# generating an image
fig = plt.figure()
plt.xlabel('time [s]')
plt.ylabel('position [m]')
plt.title('Forced Oscillator')

# set oscillator
oscill = ForcedOscillator(10, 9.5)

# simulation with different increment 
for h in {0.0001} :

    # set position
    oscill.set_position([0, 0, 0], [0, 0, 0])

    while (t < tmax) :
        # updating data for plotting
        time.append(t)
        coord.append(oscill.get_coord_x())

        # ode solver
        oscill.euler(h, t)
        t += h

    # plot position
    plt.plot(time, coord)

    t = 0
    time = []
    coord = []

print("Simulation ended!") 
print("Image generated!")
fig.savefig('ForcedOscillator.jpg', bbox_inches='tight')
plt.show()
