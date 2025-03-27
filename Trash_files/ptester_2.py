import numpy as np
from scipy.fft import fft, ifft, fftshift, fftfreq

# Define grid

length_x = 25
x = np.linspace(0, 25, 500)
dx = x[1] - x[0]

# Define function v(x) (Gaussian)
v_x = np.exp(-x**2)


dkx = length_x / (500 - 1)
kx = np.fft.fftfreq(int(500), d=dkx)

def gaussian (k):
    epsilon  = 2
    sigma  = 1
    return epsilon * sigma**3 * np.pi**(3/2) * np.exp(-((2*np.pi* k *sigma)**2/4.0))

rho = []
weight  = []
for i in range(len(kx)):
    k = kx[i] 
    weight.append(gaussian(k))
    
    rho.append(3.0)
rho =np.array(rho)
weight = np.array(weight)

print(weight)


rho_k = np.zeros(len(kx))
rho_k = fft(rho)
print(rho_k)
conv =  rho_k*weight

intgrl = ifft(conv).real


print(intgrl)
