import numpy as np
import matplotlib.pyplot as plt
import threading
import time

class Orbit(threading.Thread):
    i = 0

    def __init__(self, deltaT, init_pos, init_vel, name=None, steps=None):
        if name is None:
            name = "Orbit simulation (Generic) #{0}".format(Orbit.i)
        Orbit.i += 1

        threading.Thread.__init__(self, name=name)

        if steps is None:
            steps = int(1 / deltaT)
        self.steps = steps

        self.positions = [[init_pos[0]], [init_pos[1]]]
        self.velocities = [[init_vel[0]], [init_vel[1]]]
        self.potential_energy = [-4 * np.pi**2 / (np.linalg.norm(init_pos))]
        self.kinetic_energy = [0.5 * np.linalg.norm(init_vel)**2]

        self.deltaT = deltaT
        self.name = name

    def energy(self, pos=None, vel=None):
        if pos is None:
            pos = (self.positions[0][-1], self.positions[1][-1])
        if vel is None:
            vel = (self.velocities[0][-1], self.velocities[1][-1])
        pot = -4 * np.pi**2 / (np.linalg.norm(pos))
        kin = 0.5 * np.linalg.norm(vel)**2
        return pot, kin

    def step(self):
        pass

    def run(self):
        for i in range(self.steps):
            self.step()


class RK4Orbit(Orbit):
    rk_i = 0

    def __init__(self, deltaT, init_pos, init_vel, name=None, steps=None):
        if name is None:
            name = "Orbit simulation (Runge-Kutta 4, h = {0}) #{1}".format(deltaT, RK4Orbit.rk_i)
        RK4Orbit.rk_i += 1
        super().__init__(deltaT, init_pos, init_vel, name, steps)

    def accel(self, pos):
        r3 = np.linalg.norm(pos)**3
        accel = -4 * np.pi**2 * pos / r3
        return accel

    def step(self, deltaT=None):
        if deltaT is None:
            deltaT = self.deltaT

        i = len(self.kinetic_energy) - 1
        print("{0}: Step {1}".format(self.name, i))
        pos = np.asarray([self.positions[0][-1], self.positions[1][-1]])
        vel = np.asarray([self.velocities[0][-1], self.velocities[1][-1]])

        k1_v = self.accel(pos)
        k1_r = vel

        k2_v = self.accel(pos + k1_r * deltaT / 2)
        k2_r = vel + deltaT / 2 * k1_v

        k3_v = self.accel(pos + k2_r * deltaT / 2)
        k3_r = vel + deltaT / 2 * k2_v

        k4_v = self.accel(pos + k3_r * deltaT)
        k4_r = vel + deltaT * k3_v

        new_vel = vel + (deltaT / 6) * (k1_v + 2*(k2_v + k3_v) + k4_v)
        new_pos = pos + (deltaT / 6) * (k1_r + 2*(k2_r + k3_r) + k4_r)

        self.velocities[0].append(new_vel[0])
        self.velocities[1].append(new_vel[1])

        self.positions[0].append(new_pos[0])
        self.positions[1].append(new_pos[1])

        pot, kin = self.energy()
        self.potential_energy.append(pot)
        self.kinetic_energy.append(kin)

        return new_vel, new_pos


class ECOrbit(Orbit):
    ec_i = 0

    def __init__(self, deltaT, init_pos, init_vel, name=None, steps=None):
        if name is None:
            name = "Orbit simulation (Euler-Cromer, h = {0}) #{1}".format(deltaT, ECOrbit.ec_i)
        ECOrbit.ec_i += 1
        super().__init__(deltaT, init_pos, init_vel, name, steps)

    def step(self, deltaT=None):
        if deltaT is None:
            deltaT = self.deltaT
        
        i = len(self.positions[0]) - 1

        print("{0}: Step {1}".format(self.name, i))
        
        r = np.sqrt(self.positions[0][i]**2 + self.positions[1][i]**2)

        vx = self.velocities[0][i] - deltaT * (4 * np.pi**2 * self.positions[0][i] / (r**3))
        vy = self.velocities[1][i] - deltaT * (4 * np.pi**2 * self.positions[1][i] / (r**3))

        self.velocities[0].append(vx)
        self.velocities[1].append(vy)

        x = self.positions[0][i] + deltaT * vx
        y = self.positions[1][i] + deltaT * vy

        self.positions[0].append(x)
        self.positions[1].append(y)

        return (vx, vy), (x, y)



lin200 = np.linspace(0.00001, 0.001, 2000)
lin201 = np.linspace(0.01, 0.2, 1)

ts = time.time()

orbits200 = [RK4Orbit(i, (1, 0), (0, 2*np.pi), steps=int(1/i)*3) for i in lin200]
orbits201 = [RK4Orbit(i, (1, 0), (0, 2*np.pi), steps=int(1/i)*3) for i in lin201]

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




















