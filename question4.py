import numpy as np
import matplotlib.pyplot as plt
from orbits import RK4DualOrbit
import utilities

init_pos_e = (utilities.semimajor['earth'] * (1 - utilities.eccentricity['earth']), 0)
init_vel_e = (0, utilities.v_max('earth'))
m_e = utilities.masses_e['earth']

init_pos_m = (utilities.semimajor['mars'] * (1 - utilities.eccentricity['mars']), 0)
init_vel_m = (0, utilities.v_max('mars'))
m_m = utilities.masses_e['mars']

m_s = utilities.masses_e['sun']

def part_1(step, steps):

    earth = RK4DualOrbit(step, init_pos_e, init_vel_e, m_e, 'Earth')
    mars = RK4DualOrbit(step, init_pos_m, init_vel_m, m_m, 'Mars')

    earth.link(mars)
    mars.link(earth)

    for i in range(steps):
        earth.step()
        mars.step()

    max_x = 1.1 * np.linalg.norm(mars.positions[0], np.inf)
    max_y = 1.1 * np.linalg.norm(mars.positions[1], np.inf)

    plt.plot(mars.positions[0], mars.positions[1], 'k.', markersize=1, label="Mars")
    plt.plot(earth.positions[0], earth.positions[1], 'k.', markersize=1, label="Earth")
    plt.axes().set_aspect('equal', 'box')

    plt.xlim([-max_x, max_x])
    plt.ylim([-max_y, max_y])

    plt.axvline(x=0, linestyle='-', color='k', lw=0.4)
    plt.axhline(y=0, linestyle='-', color='k', lw=0.4)

    #plt.plot(0, 0, 'y.', markersize=30)

    sun = plt.Circle((0,0), 0.1, color='y')
    plt.axes().add_artist(sun)

    plt.legend()

    plt.show()

part_1(0.001, int(10/0.001))