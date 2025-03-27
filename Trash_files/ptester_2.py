import numpy as np
import matplotlib.pyplot as plt

# Define the Gaussian function, but only for x > 0
def gaussian(x, sigma):
    return 2* np.pi * np.exp(-x**2 / (sigma**2))

# Define analytical Fourier Transform for comparison (general case)
def gaussian_ft(f, sigma):
    return 2* np.pi * sigma * np.sqrt( np.pi) * np.exp(- (np.pi * sigma * f)**2)

# Define the domain (only positive x values)
N = 500  # Number of points
L = 25.0  # Upper limit of x
x = np.linspace(0, L, N)  # Only positive x values
sigma = 1.0  # Width of Gaussian

# Compute the Gaussian function
g = gaussian(x, sigma)

# Compute FFT
G_fft = np.fft.fft(g)
#G_fft = np.fft.fftshift(G_fft)  # Shift zero frequency to center
freqs = (np.fft.fftfreq(N, d=(L)/N))  # Adjust frequency range


print (G_fft*(2*L/N))

exit(0)

# Normalize FFT output
G_fft_magnitude = np.abs(G_fft) * (2*L/N)  # Scaling factor

# Compute analytical result (this is not exact anymore, but still useful)
G_analytical = gaussian_ft(freqs, sigma)

# Plot results
plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)
plt.plot(x, g, label="Gaussian Function (x > 0)")
plt.xlabel("x")
plt.ylabel("Amplitude")
plt.legend()
plt.title("Truncated Gaussian (Only Positive x)")

plt.subplot(1, 2, 2)
plt.plot(freqs, G_fft_magnitude, label="FFT Result", color='red', linestyle='dashed')
plt.plot(freqs, G_analytical, label="Analytical Result (Full Gaussian)", color='blue', linestyle='solid')
plt.xlabel("dfdfFrequency")
plt.ylabel("Magnitude")
plt.legend()
plt.title("Fourier Transform of Truncated Gaussian")

plt.show()

