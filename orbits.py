import threading
import numpy as np
import utilities


class Orbit(threading.Thread):
    i = 0

    def __init__(self, deltaT, init_pos, init_vel, name=None, steps=None, alpha=0, beta=2):
        if name is None:
            name = "Orbit simulation (Generic) #{0}".format(Orbit.i)
        Orbit.i += 1

        threading.Thread.__init__(self, name=name)

        if steps is None:
            steps = int(1 / deltaT)
        self.steps = steps

        self.positions = [[init_pos[0]], [init_pos[1]]]
        self.velocities = [[init_vel[0]], [init_vel[1]]]

        if self.energy() is not None:
            pot, kin = self.energy()
            self.potential_energy = [pot]
            self.kinetic_energy = [kin]

        self.deltaT = deltaT
        self.name = name
        self.alpha = alpha
        self.beta = beta

    def pos(self):
        return np.asarray([self.positions[0][-1], self.positions[1][-1]])

    def vel(self):
        return np.asarray([self.velocities[0][-1], self.velocities[1][-1]])

    def accel(self, pos):
        r = np.linalg.norm(pos)
        accel = -4 * np.pi**2 * pos / r**(self.beta + 1) * (1 + self.alpha / r**2)
        return accel

    def energy(self, pos=None, vel=None):
        if pos is None:
            pos = self.pos()
        if vel is None:
            vel = self.vel()
        pot = -4 * np.pi ** 2 / (np.linalg.norm(pos))
        kin = 0.5 * np.linalg.norm(vel) ** 2
        return pot, kin

    def step(self):
        pass

    def run(self):
        for i in range(self.steps):
            self.step()


class RK4Orbit(Orbit):
    rk_i = 0

    def __init__(self, deltaT, init_pos, init_vel, name=None, steps=None, alpha=0, beta=2):
        if name is None:
            name = "Orbit simulation (Runge-Kutta 4, h = {0}) #{1}".format(deltaT, RK4Orbit.rk_i)
        RK4Orbit.rk_i += 1
        super().__init__(deltaT, init_pos, init_vel, name, steps, alpha, beta)

    def step(self, deltaT=None):
        if deltaT is None:
            deltaT = self.deltaT

        i = len(self.kinetic_energy) - 1
        print("{0}: Step {1}".format(self.name, i))
        pos = self.pos()
        vel = self.vel()

        k1_v = self.accel(pos)
        k1_r = vel

        k2_v = self.accel(pos + k1_r * deltaT / 2)
        k2_r = vel + deltaT / 2 * k1_v

        k3_v = self.accel(pos + k2_r * deltaT / 2)
        k3_r = vel + deltaT / 2 * k2_v

        k4_v = self.accel(pos + k3_r * deltaT)
        k4_r = vel + deltaT * k3_v

        new_vel = vel + (deltaT / 6) * (k1_v + 2 * (k2_v + k3_v) + k4_v)
        new_pos = pos + (deltaT / 6) * (k1_r + 2 * (k2_r + k3_r) + k4_r)

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

    def __init__(self, deltaT, init_pos, init_vel, name=None, steps=None, alpha=0, beta=2):
        if name is None:
            name = "Orbit simulation (Euler-Cromer, h = {0}) #{1}".format(deltaT, ECOrbit.ec_i)
        ECOrbit.ec_i += 1
        super().__init__(deltaT, init_pos, init_vel, name, steps, alpha, beta)

    def step(self, deltaT=None):
        if deltaT is None:
            deltaT = self.deltaT

        i = len(self.positions[0]) - 1

        print("{0}: Step {1}".format(self.name, i))

        pos = self.pos()
        vel = self.vel()

        new_vel = vel + deltaT * self.accel(pos)
        new_pos = pos + deltaT * new_vel

        self.velocities[0].append(new_vel[0])
        self.velocities[1].append(new_vel[1])

        self.positions[0].append(new_pos[0])
        self.positions[1].append(new_pos[1])

        pot, kin = self.energy()
        self.potential_energy.append(pot)
        self.kinetic_energy.append(kin)

        return new_vel, new_pos


class RK4OrbitSI(RK4Orbit):
    rksi_i = 0

    def __init__(self, deltaT, init_pos, init_vel, name=None, steps=None, alpha=0, beta=2):
        if name is None:
            name = "Orbit simulation (Runge-Kutta 4 SI, h = {0}) #{1}".format(deltaT, RK4OrbitSI.rksi_i)
        RK4OrbitSI.rksi_i += 1
        super().__init__(deltaT, init_pos, init_vel, name, steps, alpha, beta)

    def accel(self, pos):
        r = np.linalg.norm(pos)
        accel = -utilities.G * utilities.masses_kg['earth'] * pos / r**(self.beta + 1) * (1 + self.alpha / r**2)
        return accel

    def energy(self, pos=None, vel=None):
        if pos is None:
            pos = self.pos()
        if vel is None:
            vel = self.vel()
        pot = -utilities.G * utilities.masses_kg['earth'] / np.linalg.norm(pos)
        kin = 0.5 * np.linalg.norm(vel) ** 2
        return pot, kin


class RK4DualOrbit(RK4Orbit):
    rkdual_i = 0

    def __init__(self, deltaT, init_pos, init_vel, mass, name=None, steps=None):
        if name is None:
            name = "Orbit simulation (Euler-Cromer, h = {0}) #{1}".format(deltaT, RK4DualOrbit.rkdual_i)
        RK4DualOrbit.rkdual_i += 1
        self.inter_planet: RK4DualOrbit = None
        self.mass = mass
        super().__init__(deltaT, init_pos, init_vel, name, steps, alpha=0, beta=2)

    def link(self, inter_planet: Orbit):
        self.inter_planet = inter_planet

        pot, kin = self.energy()
        self.potential_energy = [pot]
        self.kinetic_energy = [kin]

    def accel(self, pos):
        if self.inter_planet is None:
            return None

        r = np.linalg.norm(pos)
        GMs = utilities.GMs
        ratio = self.inter_planet.mass / utilities.masses_e['sun']
        inter_pos = self.inter_planet.pos()
        accel = -(GMs * pos / r**3) - (ratio*GMs * (pos - inter_pos) / r**3)
        return accel

    def energy(self, pos=None, vel=None):
        if self.inter_planet is None:
            return None

        if pos is None:
            pos = self.pos()
        if vel is None:
            vel = self.vel()

        GMs = utilities.GMs
        inter_pos = self.inter_planet.pos()
        ratio = self.inter_planet.mass / utilities.masses_e['sun']
        pot = -(GMs / np.linalg.norm(pos)) - (ratio * GMs / np.linalg.norm(pos - inter_pos))
        kin = 0.5 * np.linalg.norm(vel) ** 2
        return pot, kin

    def step(self, deltaT=None):
        if self.inter_planet is None:
            return
        super().step(deltaT)

    def run(self):
        if self.inter_planet is None:
            return
        for i in range(self.steps):
            self.step()
