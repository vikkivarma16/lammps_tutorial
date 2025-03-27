import numpy as np
import matplotlib.pyplot as plt
from scipy.special import j0

# Define Gaussian Potential
def gaussian_potential(r, V0, rc, sigma):
    return V0 * np.exp(-((r - rc)**2) / (2 * sigma**2))

# Parameters
V0 = 1.0       # Amplitude of Gaussian
rc = 5.0       # Center of Gaussian
sigma = 1.0    # Width of Gaussian
k = 2.0        # Wavenumber
r_min, r_max = 0, 10  # Integration limits

# Generate r values
r_values = np.linspace(r_min, r_max, 500)
V_values = gaussian_potential(r_values, V0, rc, sigma)
j0_values = j0(k * r_values)
weighted_values = 4 * np.pi * r_values**2 * V_values * j0_values

# Plot the interaction potential
plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
plt.plot(r_values, V_values, label="Gaussian Potential V(r)")
plt.xlabel("r")
plt.ylabel("V(r)")
plt.legend()
plt.title("Positive Gaussian Potential")

# Plot the weighted function V(r) * j0(kr)
plt.subplot(1, 2, 2)
plt.plot(r_values, weighted_values, label="Weighted Function V(r) j0(kr)")
plt.axhline(0, color='k', linestyle='--')
plt.xlabel("r")
plt.ylabel("V(r) j0(kr)")
plt.legend()
plt.title("Weighted Function V(r) j0(kr)")

plt.show()

