import numpy as np
import matplotlib.pyplot as plt

def load_data(filename):
    """Loads data from a text file assuming two columns: x and y."""
    data = np.loadtxt(filename)
    x = data[:, 0]
    y = data[:, 1]
    return x, y

def plot_data(x, y, output_filename='data_plot.png'):
    """Plots the data and saves it as a high-resolution PNG file."""
    plt.figure(figsize=(8, 5))
    plt.plot(x, y, marker='o', linestyle='-', color='b')
    plt.xlabel('rho_p')
    plt.ylabel('Fraction of adsorbed polymers')
    plt.title('Polymer density vs absorptivity')
    plt.grid()
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    input_filename = 'data.txt'  # Change this if your filename is different
    x, y = load_data(input_filename)
    plot_data(x, y)
