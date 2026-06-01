import numpy as np
import matplotlib.pyplot as plt

def n_o(lmbd_um):
    return np.sqrt(
        1
        + 2.6734 / (1 - 0.01764 / lmbd_um**2)
        + 1.2290 / (1 - 0.05914 / lmbd_um**2)
        + 12.614 / (1 - 474.60 / lmbd_um**2)
    )
def n_e(lmbd_um):
    return np.sqrt(
        1
        + 2.9804 / (1 - 0.02047 / lmbd_um**2)
        + 0.5981 / (1 - 0.0666 / lmbd_um**2)
        + 8.9543 / (1 - 416.08 / lmbd_um**2)
    )
def rhs(z, Y, k0, no, ne):
    Ex, Ey = Y

    return np.array([
        1j * k0 * no * Ex,
        1j * k0 * ne * Ey
    ], dtype=complex)
def rk4_step(f, z, Y, h, *args):

    k1 = f(z, Y, *args)
    k2 = f(z + h / 2, Y + h * k1 / 2, *args)
    k3 = f(z + h / 2, Y + h * k2 / 2, *args)
    k4 = f(z + h, Y + h * k3, *args)

    return Y + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6

lambda_nm = 550.0
lambda_um = lambda_nm / 1000.0
lambda_m = lambda_nm * 1e-9

d_um = 100.0
d_m = d_um * 1e-6

theta_deg = 45
theta = np.radians(theta_deg)

no = n_o(lambda_um)
ne = n_e(lambda_um)

delta_n = ne - no

print(f"λ = {lambda_nm} nm")
print(f"n_o = {no:.6f}")
print(f"n_e = {ne:.6f}")
print(f"Δn = {delta_n:.6f}")

k0 = 2 * np.pi / lambda_m

N = 2000

z = np.linspace(0, d_m, N)
h = z[1] - z[0]

Y = np.zeros((N, 2), dtype=complex)

Y[0] = [
    np.cos(theta),
    np.sin(theta)
]

for i in range(N - 1):
    Y[i + 1] = rk4_step(
        rhs,
        z[i],
        Y[i],
        h,
        k0,
        no,
        ne
    )

Ex_num = Y[-1, 0]
Ey_num = Y[-1, 1]

Ex_exact = np.cos(theta) * np.exp(1j * k0 * no * d_m)
Ey_exact = np.sin(theta) * np.exp(1j * k0 * ne * d_m)

error = np.sqrt(
    np.abs(Ex_num - Ex_exact)**2 +
    np.abs(Ey_num - Ey_exact)**2
)

print(f"RK4 error = {error:.3e}")
delta_phi = k0 * (ne - no) * d_m

print(f"Phase delay Δφ = {delta_phi:.4f} rad")
print(f"Δφ / π = {delta_phi/np.pi:.4f}")
t = np.linspace(0, 2*np.pi, 1000)

x_in = np.cos(theta) * np.cos(t)
y_in = np.sin(theta) * np.cos(t)

x_out = np.real(Ex_num * np.exp(1j * t))
y_out = np.real(Ey_num * np.exp(1j * t)

lambda_range_nm = np.linspace(450, 650, 400)
lambda_range_um = lambda_range_nm / 1000

no_arr = n_o(lambda_range_um)
ne_arr = n_e(lambda_range_um)

delta_n_arr = ne_arr - no_arr

delta_phi_arr = (
    2*np.pi /
    (lambda_range_nm * 1e-9)
    * delta_n_arr
    * d_m
)

plt.figure(figsize=(8,5))
plt.plot(lambda_range_nm, no_arr, label='n_o')
plt.plot(lambda_range_nm, ne_arr, label='n_e')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Refractive index')
plt.title('LiNbO3 refractive indices')
plt.grid(True)
plt.legend()

plt.figure(figsize=(8,5))
plt.plot(lambda_range_nm, delta_n_arr)
plt.xlabel('Wavelength (nm)')
plt.ylabel('Δn')
plt.title('Birefringence')
plt.grid(True)

plt.figure(figsize=(8,5))
plt.plot(lambda_range_nm, delta_phi_arr)
plt.xlabel('Wavelength (nm)')
plt.ylabel('Δφ (rad)')
plt.title('Phase retardation')
plt.grid(True)

plt.figure(figsize=(12,5))

plt.subplot(1,2,1)
plt.plot(x_in, y_in)
plt.axis('equal')
plt.grid(True)
plt.xlabel('Ex')
plt.ylabel('Ey')
plt.title('Input polarization')

plt.subplot(1,2,2)
plt.plot(x_out, y_out)
plt.axis('equal')
plt.grid(True)
plt.xlabel('Ex')
plt.ylabel('Ey')
plt.title('Output polarization')

plt.tight_layout()
plt.show()
