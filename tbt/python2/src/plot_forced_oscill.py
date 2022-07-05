
# author:          Lorenzo Liuzzo
# email:           lorenzoliuzzo@outlook.com
# description:     Plot of a forced oscillator using ODE's numerical methods.
# last updated:    02/07/2022


import matplotlib.pyplot as plt
from oscillators import ForcedOscillator
from position import Position


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
oscill = ForcedOscillator(10, 5)

# simulation with different increment 
for h in {0.0001} :

    # set position
    pos = Position([0, 0, 0], [0, 0, 0])

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
fig.savefig('ForcedOscillator.jpg', bbox_inches='tight')
plt.show()
