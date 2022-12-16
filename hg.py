import numpy as np
from scipy.optimize import fsolve as fsolve
import matplotlib.pyplot as plt


# Location Parameters
diameter = 0.1016						# m
isNozzle = False						# Is the region you're interested in in the nozzle? This would be anything not in the straight portion of the chamber.
T_wg_Tc_ratio = 0.4

# Other parameters
area = np.pi * diameter**2 / 4    		# m^2

# Chamber Parameters
m_dot = 0.958                           # kg/s
P_chamber = 500 * 6894.76				# Pa
T_chamber = 3030						# K

# Geometry
d_throat = 0.02353                      # m
A_throat = np.pi * d_throat**2 / 4      # m^2
nozzle_radius_of_curvature = 0.0183261	# m

# Fluid Parameters
g = 1.24
M = 23.9								# kg / kmol
density = 3.2662						# kg / m^3
mu = 8.686258984659215e-5				# kg / m·s
Cp = 1778.3692185308112					# J/(kg·K)
k = 0.2613354950693967					# W/(m·K)

Pr = mu * Cp / k
C_star = P_chamber * A_throat / m_dot	# m/s

# Ratio duct area to critical duct area (location where M=1) for isentropic flow. M is the mach number.
def isentropic_A_A_star(M):
    return ((2 / (g+1)) * (1 + 0.5 * (g-1) * M**2)) ** ((g+1)/2/(g-1)) / M

if (isNozzle):
	func = lambda M : (diameter/d_throat)**2 - isentropic_A_A_star(M)
	mach = fsolve(func, 0.5)[0]
	sigma = (0.5 * T_wg_Tc_ratio * (1 + mach**2 * (g-1)/2) + 0.5)**(-0.68) * (1 + mach**2 * (g-1)/2)**(-0.12)
else:
	sigma = 1

hg = (0.026 / d_throat**0.2) * (mu**0.2 * Cp / Pr**0.6) * (P_chamber / C_star)**0.8 * (d_throat / nozzle_radius_of_curvature)**0.1 * (A_throat / area)**0.9 * sigma
print(f"Sigma: {sigma}")
print(f"h_g: {hg}")


# The follow parameters are from Sample Calculation 4-3 in Huzel and Huang

# d_throat = 0.63246                      # m
# nozzle_radius_of_curvature = 0.297434	# m

# P_chamber = 1000 * 6894.76				# Pa
# T_chmaber = 3590						# K

# g = 1.222
# M = 22.5
# mu = 0.000074646301224					# kg / m·s
# Cp = 2030								# J/(kg·K)

# Pr = 0.816

# C_star = 1725.168
# sigma = 1.05