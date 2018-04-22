import numpy as np

G = 6.67408E-11

# GMs in units AU and yr
GMs = 4 * np.pi ** 2

# Masses in units of the mass of the Earth
masses_e = {'sun': 333480,
            'mercury': 0.0553,
            'venus': 0.8150,
            'earth': 1.000,
            'mars': 0.1074,
            'jupiter': 317.89,
            'uranus': 14.56,
            'neptune': 17.15,
            'pluto': 0.002}

masses_kg = {}
for k,v in masses_e.items():
    masses_kg[k] = v * 6.0E24

for v in masses_kg.values():
    print(v)

eccentricity = {'mercury': 0.2056,
                'venus': 0.0068,
                'earth': 0.0167,
                'mars': 0.0934,
                'jupiter': 0.0483,
                'saturn': 0.0560,
                'uranus': 0.0461,
                'neptune': 0.0100,
                'pluto': 0.2484}

# Semimajor axis in AU
semimajor = {'mercury': 0.3871,
             'venus': 0.7233,
             'earth': 1.0000,
             'mars': 1.5237,
             'jupiter': 5.2028,
             'saturn': 9.5388,
             'uranus': 19.191,
             'neptune': 30.061,
             'pluto': 39.529}

# Period in yrs
period = {'mercury': 0.2408,
          'venus': 0.6152,
          'earth': 1.0000,
          'mars': 1.8809,
          'saturn': 29.456,
          'uranus': 84.07,
          'neptune': 164.81,
          'pluto': 248.53}

radius = {'earth': 6.3781E6}

def v_max(planet):
    e = eccentricity[planet]
    a = semimajor[planet]
    Mp = masses_e[planet]
    Ms = masses_e['sun']
    ratio = Mp / Ms
    return np.sqrt(GMs * (1 + e) / (a * (1 - e)) * (1 + ratio))

def v_max_SI(M, m, Rp, Ra):
    e = (Ra - Rp) / (Ra + Rp)
    ratio = m / M
    a = Rp / (1 - e)
    return np.sqrt(G * M * (1 + e) / (a * (1 - e)) * (1 + ratio))