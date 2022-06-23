from matplotlib import pyplot as plt
from ode import ForcedOscillator

# oscillator setting
oscill = ForcedOscillator(10, 5)
oscill.set_position([0, 0, 0], [0, 0, 0])

t = 0
# set tmax and the increment of t for each cicle
tmax = 100
h = 0.001

time = []
x = []
vel = []

# simulation
oscill.print_position()
while (t < tmax) :
    # updating data for plotting
    time.append(t)
    x.append(oscill.get_coord_x())
    vel.append(oscill.get_vel_x())
    
    # ode solver
    oscill.simpleEuler(h)
    t += h
    
    #animation 
    plt.plot(time, x)
    plt.pause(0.01)
plt.show()
print("Simulazione finita!")

# generating an image
fig = plt.figure(figsize = (15, 15))
plt.title('Forced Oscillator', fontsize = 15)

# plot position
sub1 = plt.subplot(2, 1, 1)
sub1.plot(time, x)
sub1.set_xlabel('time [s]')
sub1.set_ylabel('position [m]')
sub1.set_title('Position')

# plot velocity
sub2 = plt.subplot(2, 1, 2)
sub2.plot(time, vel, 'g')
sub2.set_xlabel('time [s]')
sub2.set_ylabel('velocity [m/s]')
sub2.set_title('Velocity')

fig.savefig('ForcedOscillator.jpg', bbox_inches='tight', dpi=150)
print("Immagine creata!")
plt.show()
