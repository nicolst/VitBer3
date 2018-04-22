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
    mercury = RK4Orbit(step_size, start_pos, start_vel, steps=steps)
    mercury.start()
    mercury.join()

    max_x = 1.1 * np.linalg.norm(mercury.positions[0], np.inf)
    max_y = 1.1 * np.linalg.norm(mercury.positions[1], np.inf)

    plt.figure(1)

    plt.plot(mercury.positions[0], mercury.positions[1], 'k.', markersize=0.4)
    plt.axes().set_aspect('equal', 'box')
    plt.xlim([-max_x, max_x])
    plt.ylim([-max_y, max_y])

    plt.axvline(x=0, linestyle='-', color='k', lw=0.4)
    plt.axhline(y=0, linestyle='-', color='k', lw=0.4)

    plt.plot(0, 0, 'y.', markersize=30)

    plt.xlabel(r"$x$ (AU)")
    plt.ylabel(r"$y$ (AU)")

    plt.title("Orbit of Mercury")

    plt.tight_layout()

    plt.figure(2)
    plt.title("Velocity of Mercury")
    plt.plot(np.arange(steps)*step_size, [np.sqrt(mercury.velocities[0][i]**2 + mercury.velocities[1][i]**2) for i in range(steps)], 'k-')

    plt.ylabel(r"$|v|$ (AU/yr)")
    plt.xlabel(r"$t$ (yr)")

    plt.tight_layout()

    plt.figure(3)
    total_energy = np.asarray(mercury.kinetic_energy) + np.asarray(mercury.potential_energy)
    plt.plot(np.arange(steps)*step_size, total_energy[:-1], 'k-', label="Total energy")
    plt.plot(np.arange(steps)*step_size, mercury.kinetic_energy[:-1], 'k--', label="Kinetic energy")
    plt.plot(np.arange(steps)*step_size, mercury.potential_energy[:-1], 'k:', label="Potential energy")

    plt.title("Total, kinetic and potential energy for orbit of Mercury")
    plt.ylabel(r"$E/n$ (AU$^2$ / yr$^2$)")
    plt.xlabel(r"$t$ (yr)")

    plt.tight_layout()
    plt.legend()
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
part_1(step_size, int(3 * period / step_size))
#part_2(step_size, int(0.2 / step_size), beta=2.5)













