import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

k = 398    # W/m*K
rho = 8960  # kg/m^3
c = 385     # J/kg*K
a = k / (c * rho)

t_final = 5        # second
x_final = 0.071     # height of housing in meters
delta_t = 5e-6      # Time step in seconds, NOTE: Delta_t must be of the same order of magnitude as delta_x squared, preferably even smaller.
num_x_points = 100  # Number of subdivisions
T_combustion = 3029	# Kelvin
T_external = 300	# Kelvin
hg = 720			# W / (m^2 · K)
area = 4.03225e-5	# m^2

# --------------------------------- [Computation] --------------------------------- #

num_t_points = round(t_final / delta_t)
delta_x = x_final / num_x_points

T = np.zeros((num_t_points, num_x_points + 1))	# The + 1 is for the fake boundary condition point on the outer wall
T[:, 0:-1] = 300
T[:, -1] = T_external

for i in range(1, num_t_points):
	# dT/dt = a * d^2T/dx^2
	second_dT_dx = (T[i-1, 2:] - 2 * T[i-1, 1:-1] + T[i-1, :-2]) / (delta_x ** 2)
	dT_dt = a * second_dT_dx
	T[i, 1:-1] = T[i-1, 1:-1] + dT_dt * delta_t

	# Compute the dT_dt for the hot-wall side seperately by applying Newton's Law of Cooling and Fourier's Law of Heat Conduction
	# dT/dt = A*hg/C * (Tc - T) - A*k/C/delta_x * (T - T_1)
	heat_capacity = c * (rho * area * delta_x)
	dT_dt = (area * hg / heat_capacity) * (T_combustion - T[i-1, 0]) - (area * k / heat_capacity) * (T[i-1, 0] - T[i-1, 1]) / delta_x
	T[i, 0] = T[i-1, 0] + dT_dt * delta_t

T = T[:, 0:-1]	# Remove the fake boundary condition point

# --------------------------------- [Plotting] --------------------------------- #
fig, ax = plt.subplots()
line, = ax.plot(np.linspace(0, x_final, num_x_points) * 1000, T[0] - 273)
ax.set_ylim(0, 2000)
ax.set_xlabel("Location (mm)")
ax.set_ylabel("Temperature (°C)")

# Init only required for blitting to give a clean slate.
def init():
    line.set_ydata([])
    return line,

def animate(i):
	line.set_ydata(T[i] - 273)
	return line,

ani = animation.FuncAnimation(fig, animate, init_func=init, frames=range(0, num_t_points, 1000), interval=1, blit=True, repeat=False)

# writer = animation.PillowWriter(fps=60)
# dir = "C:\\Users\\Jiyansh User\\Desktop"
# ani.save(f"{dir}\\1_second_run_copper.gif", writer=writer)
plt.show()