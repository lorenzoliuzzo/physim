
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

print("This is a simulation of a forced oscillator.\n")

# set oscillator
oscill = ForcedOscillator(10, 5)

# generating an image
fig = plt.figure(figsize = (10, 10))
plt.title('Forced Oscillator', fontsize = 15)

# simulation with different increment 
for h in {0.01, 0.001, 0.0001} :

    # set position
    oscill.set_position([0, 0, 0], [0, 0, 0])
    while (t < tmax) :
        # updating data for plotting
        time.append(t)
        coord.append(oscill.get_coord_x())

        # ode solver
        oscill.rk4(h, t)
        t += h

    t = 0

    # plot position
    plt.xlabel('time [s]')
    plt.ylabel('position [m]')
    plt.title('Position')
    plt.plot(time, coord)

    time = []
    coord = []
    
print("Simulazione finita!")
print("Immagine creata!")
fig.savefig('ForcedOscillator.jpg', bbox_inches='tight', dpi=150)
plt.show()
