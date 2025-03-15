import numpy as np
import matplotlib.pyplot as plt

def parse_lammps_chunk_data(filename):
    time_data = {}
    with open(filename, 'r') as file:
        lines = file.readlines()
        timestep = None
        i = 0
        threshold = 0
        for line in lines:
            parts = line.strip().split()
            if len(parts) == 2 and i == 0:
                timestep = int(parts[0]) * 0.005
                time_data[timestep] = []
                threshold = int(parts[1])
                i += 1
            elif len(parts) == 2 and i <= threshold:
                chunk, count = map(int, parts)
                time_data[timestep].append((chunk, count))
                i += 1
                if i > threshold:
                    i = 0
    return time_data

def calculate_adsorbed_fraction(time_data, bin_volume, dx, sticky_particles_bin_range):
    adsorbed_fraction = {}
    total_adsorbed = 0
    for timestep, data in time_data.items():
        bin_sums = {}
        bin_counts = {}
        for bin_id, count in data:
            if bin_id not in bin_sums:
                bin_sums[bin_id] = 0
                bin_counts[bin_id] = 0
            bin_sums[bin_id] += count
            bin_counts[bin_id] += 1
        
        bin_positions = []
        densities = []
        for bin_id in sorted(bin_sums.keys()):
            avg_count = bin_sums[bin_id] / bin_counts[bin_id]
            density = avg_count / bin_volume
            bin_positions.append(bin_id * dx)
            densities.append(density)
        
        densities = np.array(densities)
        total_particles = np.sum(densities)
        central_bin = int(0.5 * len(bin_positions))
        sticky_particles = sum(densities[i + central_bin] + densities[-i + central_bin] for i in range(sticky_particles_bin_range))
        fraction = sticky_particles / total_particles if total_particles else 0
        adsorbed_fraction[timestep] = fraction
        total_adsorbed += fraction
    
    time_averaged_adsorptivity = total_adsorbed / len(adsorbed_fraction) if adsorbed_fraction else 0
    print(f"Time-averaged adsorptivity: {time_averaged_adsorptivity}")
    return adsorbed_fraction

def save_adsorbed_fraction_plot(adsorbed_fraction, filename='adsorbed_fraction.png'):
    times = sorted(adsorbed_fraction.keys())
    fractions = [adsorbed_fraction[t] for t in times]
    
    plt.figure(figsize=(8, 5))
    plt.plot(times, fractions, marker='o', linestyle='-', color='b')
    plt.xlabel('Time')
    plt.ylabel('Fraction of Adsorbed Particles')
    plt.title('Adsorbed Fraction Over Time')
    plt.grid()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    filename = 'bin_particles.lammpstrj'
    xmin, xmax = -6.65, 6.65
    ymin, ymax = -10, 10
    zmin, zmax = -10, 10
    bin_division = 0.05
    total_volume = (xmax - xmin) * (ymax - ymin) * (zmax - zmin)
    relative_bin_size = bin_division * (xmax - xmin)
    dx = relative_bin_size
    bin_volume = bin_division * total_volume
    sticky_particles_bin_range = 3
    
    time_data = parse_lammps_chunk_data(filename)
    adsorbed_fraction = calculate_adsorbed_fraction(time_data, bin_volume, dx, sticky_particles_bin_range)
    save_adsorbed_fraction_plot(adsorbed_fraction)
