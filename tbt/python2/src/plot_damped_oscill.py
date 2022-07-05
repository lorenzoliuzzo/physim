
# author:          Lorenzo Liuzzo
# email:           lorenzoliuzzo@outlook.com
# description:     Plot of an Damped oscillator using ODE's numerical methods.
# last updated:    02/07/2022

import matplotlib.pyplot as plt
from oscillators import DampedOscillator
from position import Position


t = 0
tmax = 100

time = []
coord = []

print("This is a simulation of a damped oscillator.")

# generating an image
fig = plt.figure()
plt.xlabel('time [s]')
plt.ylabel('position [m]')
plt.title('Damped Oscillator')

# set oscillator
oscill = DampedOscillator(3., 1 / 30)

# simulation with different increment 
for h in {0.0001} :

    # set position
    pos = Position([3., 0, 0], [0, 0, 0])

    while (t < tmax) :
        # updating data for plotting
        time.append(t)
        coord.append(pos.get_coord_x())

        # ode solver
        pos.set_position(oscill.rk4(pos, h))
        t += h

    # plot position
    plt.plot(time, coord)

    t = 0
    time = []
    coord = []

print("Simulation ended!") 
print("Image generated!")
fig.savefig('DampedOscillator.jpg', bbox_inches='tight')
plt.show()
