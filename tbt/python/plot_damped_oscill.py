
# author:          Lorenzo Liuzzo
# email:           lorenzoliuzzo@outlook.com
# description:     Plot of an Damped oscillator using ODE's numerical methods.
# last updated:    24/06/2022

import matplotlib.pyplot as plt
from oscillators import DampedOscillator


t = 0
tmax = 11

time = []
coord = []

print("This is a simulation of a damped oscillator.")

# generating an image
fig = plt.figure()
plt.xlabel('time [s]')
plt.ylabel('position [m]')
plt.title('Damped Oscillator')

# set oscillator
oscill = DampedOscillator(3.1623, 1 / 10)

# simulation with different increment 
for h in {0.0001} :

    # set position
    oscill.set_position([0.196, 0, 0], [0, 0, 0])

    while (t < tmax) :
        # updating data for plotting
        time.append(t)
        coord.append(oscill.get_coord_x())

        # ode solver
        oscill.euler(h)
        t += h

    # plot position
    plt.plot(time, coord)

    t = 0
    time = []
    coord = []

print("Simulation ended!") 
print("Image generated!")
fig.savefig('es_pol.jpg', bbox_inches='tight')
plt.show()
