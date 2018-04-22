import numpy as np
import matplotlib.pyplot as plt
from orbits import RK4Orbit
import utilities

vmax = utilities.v_max('mercury')
period = utilities.period['mercury']
perihelion = utilities.semimajor['mercury'] * (1 - utilities.eccentricity['mercury'])
start_pos = (perihelion, 0)
start_vel = (0, vmax)


# steps=int(1 / step_size * period) + 1 for 1 yr
def part_1(step_size, steps):
    test = RK4Orbit(step_size, start_pos, start_vel, steps=steps)
    test.start()
    test.join()

    max_x = 1.1 * np.linalg.norm(test.positions[0], np.inf)
    max_y = 1.1 * np.linalg.norm(test.positions[1], np.inf)

    plt.plot(test.positions[0], test.positions[1], 'k.', markersize=2)
    plt.axes().set_aspect('equal', 'box')
    plt.xlim([-max_x, max_x])
    plt.ylim([-max_y, max_y])

    plt.axvline(x=0, linestyle='-', color='k', lw=0.4)
    plt.axhline(y=0, linestyle='-', color='k', lw=0.4)

    plt.plot(0, 0, 'y.', markersize=30)

    plt.show()


def part_2(step_size, steps, alpha=0, beta=2):
    test = RK4Orbit(step_size, start_pos, start_vel, steps=steps, alpha=alpha, beta=beta)
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

    plt.plot(0, 0, 'y.', markersize=30)

    plt.show()

step_size = 0.0001
#part_1(step_size, int(1 / step_size * period) + 1)
part_2(step_size, int(0.2 / step_size), beta=2.5)













