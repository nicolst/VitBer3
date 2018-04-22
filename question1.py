import numpy as np
import matplotlib.pyplot as plt
import time
from orbits import ECOrbit, RK4Orbit

lin200 = np.linspace(0.00001, 0.001, 2000)
lin201 = np.linspace(0.01, 0.2, 1)

ts = time.time()

orbits200 = [RK4Orbit(i, (1, 0), (0, 2 * np.pi), steps=int(1 / i) * 3) for i in lin200]
orbits201 = [RK4Orbit(i, (1, 0), (0, 2 * np.pi), steps=int(1 / i) * 3) for i in lin201]

for orbit in orbits200:
    orbit.start()

for orbit in orbits201:
    orbit.start()

for orbit in orbits200:
    orbit.join()

for orbit in orbits201:
    orbit.join()

tf = time.time()
print("Time taken: {0}".format(tf - ts))

plt.figure(1)
plt.title("RK4, linspace 201")
plt.plot(lin201, [(orbit.potential_energy[-1] - orbit.potential_energy[0]
                   + orbit.kinetic_energy[-1] - orbit.kinetic_energy[0]) for orbit in orbits201])

plt.figure(2)
plt.title("RK4, linspace 200")
plt.plot(lin200, [(orbit.potential_energy[-1] - orbit.potential_energy[0]
                   + orbit.kinetic_energy[-1] - orbit.kinetic_energy[0]) for orbit in orbits200])

plt.show()

"""
test = RK4Orbit(1/200, (1, 0), (0, 2*np.pi))
test2 = ECOrbit(1/10000, (1, 0), (0, 2*np.pi))

plt.plot(test.positions[0], test.positions[1], 'k.', markersize=2)
lims = (-1, 1)

plt.xticks([-1, -0.5, 0, 0.5, 1])
plt.yticks([-1, -0.5, 0, 0.5, 1])

plt.ylabel(r"$y$ (AU)")
plt.xlabel(r"$x$ (AU)")

plt.title(test.name)

plt.axhline(y=0, color='k', marker='--')
plt.axvline(x=0, color='k', marker='--')

#plt.xlim(lims)
#plt.ylim(lims)
plt.axes().set_aspect('equal', 'box')
plt.show()
"""

"""
data = []
plt.ion()
for i in range(20):
    test.step()
fig, = plt.plot(test.positions[0][-20:], test.positions[1][-20:], 'k.', markersize=2)
plt.xticks([-1, -0.5, 0, 0.5, 1])
plt.yticks([-1, -0.5, 0, 0.5, 1])

plt.ylabel(r"$y$ (AU)")
plt.xlabel(r"$x$ (AU)")

plt.title(test.name)
plt.axes().set_aspect('equal', 'box')

plt.draw()
plt.pause(0.01)
i = 0
while True:
    test.step()
    fig.set_xdata(test.positions[0][-60:])
    fig.set_ydata(test.positions[1][-60:])
    plt.xticks([-1, -0.5, 0, 0.5, 1])
    plt.yticks([-1, -0.5, 0, 0.5, 1])

    plt.ylabel(r"$y$ (AU)")
    plt.xlabel(r"$x$ (AU)")

    plt.title(test.name)
    plt.axes().set_aspect('equal', 'box')
    plt.draw()
    plt.pause(0.01)

    i += 1


plt.show(block=True)
"""

"""
rks = []
ecs = []
for i in range(10):
    rks.append(RK4Orbit(0.0001, (1, 0), (0, (i/3)*np.pi)))
    ecs.append(ECOrbit(0.0001, (1, 0), (0, (i/3)*np.pi)))

for i in range(10):
    rks[i].start()
    ecs[i].start()

for i in range(10):
    rks[i].join()
    ecs[i].join()

for i in range(10):
    plt.figure(1)
    plt.plot(rks[i].positions[0], rks[i].positions[1], '.', label=rks[i].name)
    plt.figure(2)
    plt.plot(ecs[i].positions[0], ecs[i].positions[1], '.', label=ecs[i].name)

plt.figure(1)
plt.title("Runge-Kutta")
plt.legend()

plt.figure(2)
plt.title("Euler-Cromer")
plt.legend()

plt.show()
"""
