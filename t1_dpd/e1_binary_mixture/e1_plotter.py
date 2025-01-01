import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Function to read x and Δa from a file
def read_x_vs_da(filename):
    data = np.loadtxt(filename, comments='#')
    x = data[:, 0]  # First column is x
    delta_a = data[:, 1]  # Second column is Δa
    return x, delta_a

# Function to calculate χ
def calculate_chi(x):
    # Avoid division by zero and invalid log arguments
    valid_mask = (x > 0) & (x < 1) & (1 - 2 * x != 0)
    chi = np.zeros_like(x)
    chi[valid_mask] = np.log((1 - x[valid_mask]) / x[valid_mask]) / (1 - 2 * x[valid_mask])
    return chi[valid_mask], valid_mask

# Linear function for curve fitting
def linear_func(x, m, c):
    return m * x + c

# Function to plot χ vs Δa with linear fit
def plot_chi_vs_delta_a_with_fit(delta_a, chi, output_file='chi_vs_delta_a_fit.png'):
    plt.figure(figsize=(10, 6))
    
    # Scatter plot of data points
    plt.scatter(delta_a, chi, color='b', label='Data Points')
    
    # Perform linear curve fitting
    params, _ = curve_fit(linear_func, delta_a, chi)
    m, c = params
    
    # Generate fitted line
    delta_a_fit = np.linspace(min(delta_a), max(delta_a), 500)
    chi_fit = linear_func(delta_a_fit, m, c)
    
    # Plot the fitted line
    plt.plot(delta_a_fit, chi_fit, color='r', linestyle='--', label=f'Linear Fit: $\\chi = {m:.3f} \\Delta a + {c:.3f}$')
    
    # Labels and Title
    plt.xlabel(r'$\Delta a$')
    plt.ylabel(r'$\chi$')
    plt.title(r'$\chi$ vs $\Delta a$ with Linear Fit')
    plt.legend()
    plt.grid(True)
    
    # Save the plot
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot with linear fit saved as {output_file}")

# Main execution
if __name__ == '__main__':
    filename = 'x_vs_delta_a.txt'  # Replace with your file path
    
    # Read data
    x, delta_a = read_x_vs_da(filename)
    
    # Calculate χ
    chi, valid_mask = calculate_chi(x)
    delta_a_valid = delta_a[valid_mask]
    
    # Plot and save the graph with linear fit
    plot_chi_vs_delta_a_with_fit(delta_a_valid, chi)

