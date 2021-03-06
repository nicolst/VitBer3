import numpy as np
import matplotlib.pyplot as plt
from orbits import RK4OrbitSI
import utilities

r_earth = 6.3781E6
Ms = 720

init_pos = np.asarray([322.0E3 + r_earth, 0])
init_vel = np.asarray([0, np.sqrt(utilities.G * utilities.masses_kg['earth'] / init_pos[0])])


def question_1(step):
    test = RK4OrbitSI(step, init_pos, init_vel, steps=int((90 * 60) / step))
    test.start()
    test.join()

    max_x = 1.1 * np.linalg.norm(test.positions[0], np.inf)
    max_y = 1.1 * np.linalg.norm(test.positions[1], np.inf)

    plt.plot(test.positions[0], test.positions[1], 'k.', markersize=1)
    plt.axes().set_aspect('equal', 'box')
    plt.xlim([-max_x, max_x])
    plt.ylim([-max_y, max_y])

    plt.axvline(x=0, linestyle='-', color='k', lw=0.4)
    plt.axhline(y=0, linestyle='-', color='k', lw=0.4)

    earth = plt.Circle((0, 0), r_earth, color='b')
    plt.axes().add_artist(earth)

    plt.title("Low Earth Orbit")

    plt.xlabel(r"$x$ / (m)")
    plt.ylabel(r"$y$ / (m)")

    plt.gca().get_yaxis().get_major_formatter().set_powerlimits((0, 0))
    plt.gca().get_xaxis().get_major_formatter().set_powerlimits((0, 0))

    plt.show()


def question_2(step):
    wanted_apoapsis = 35680E3 + r_earth
    wanted_periapsis = init_pos[0]
    wanted_eccentricity = (wanted_apoapsis - wanted_periapsis) / (wanted_apoapsis + wanted_periapsis)
    max_v = utilities.v_max_SI(utilities.masses_kg['earth'], Ms, wanted_periapsis, wanted_apoapsis)
    semi_major = wanted_periapsis / (1 - wanted_eccentricity)
    period = 2 * np.pi * np.sqrt(semi_major**3 / (utilities.G * utilities.masses_kg['earth']))

    # Approximately one quarter orbit in LEO
    test = RK4OrbitSI(step, init_pos, init_vel, steps=int((90 * 60) / step))
    test.run()

    # Execute first burn
    first_burn_pos = test.pos()
    prev_vel = test.vel()
    vel_ratio = max_v / np.linalg.norm(prev_vel)
    test.velocities[0][-1] *= vel_ratio
    test.velocities[1][-1] *= vel_ratio

    test.steps = int((period / 2) / step)
    test.run()

    # Execute second burn
    second_burn_pos = test.pos()
    prev_vel = test.vel()
    r = np.linalg.norm(second_burn_pos)
    new_vel = np.sqrt(utilities.G * utilities.masses_kg['earth'] / r)
    vel_ratio = new_vel / np.linalg.norm(prev_vel)
    test.velocities[0][-1] *= vel_ratio
    test.velocities[1][-1] *= vel_ratio

    period = 2 * np.pi * r / new_vel
    test.steps = int(period / step)
    test.run()
    print("got here")
    max_x = 1.1 * np.linalg.norm(test.positions[0], np.inf)
    max_y = 1.1 * np.linalg.norm(test.positions[1], np.inf)

    plt.plot(test.positions[0], test.positions[1], 'k.', markersize=0.5)
    plt.axes().set_aspect('equal', 'box')
    plt.xlim([-max_x, max_x])
    plt.ylim([-max_y, max_y])

    plt.axvline(x=0, linestyle='-', color='k', lw=0.4)
    plt.axhline(y=0, linestyle='-', color='k', lw=0.4)

    earth = plt.Circle((0, 0), r_earth, color='#0000FF77')
    plt.axes().add_artist(earth)

    plt.plot(first_burn_pos[0], first_burn_pos[1], 'ro', label="First burn")
    plt.plot(second_burn_pos[0], second_burn_pos[1], 'go', label="Second burn")

    plt.title("Hohmann transfer from LEO to GEO")
    plt.xlabel(r"$x$ / (m)")
    plt.ylabel(r"$y$ / (m)")

    plt.legend(loc="upper right")

    plt.show()

question_2(10)















