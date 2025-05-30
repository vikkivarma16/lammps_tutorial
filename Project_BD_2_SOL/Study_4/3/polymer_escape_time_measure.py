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
                timestep = int(parts[0]) * 0.005  # Time scaling factor
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

def count_escaped_particles(time_data, dx):
    escape_x_min, escape_x_max = 31, 54  # Escape region on the left
    escape_x_min_2, escape_x_max_2 = 31, 54  # Escape region on the right
    trapped_x_min, trapped_x_max = 0, 31  # Trapped region
    
    escape_fraction = {}
    total_particles_initial = None
    escape_time = None
    
    for timestep, data in time_data.items():
        escaped_particles = 0
        trapped_particles = 0
        total_particles = 0
        
        for bin_id, count in data:
            x_position = trapped_x_min + bin_id * dx  # Convert bin index to x position
            total_particles += count
            
            if escape_x_min <= x_position <= escape_x_max or escape_x_min_2 <= x_position <= escape_x_max_2:
                escaped_particles += count
            elif trapped_x_min <= x_position <= trapped_x_max:
                trapped_particles += count
        
        if total_particles_initial is None:
            total_particles_initial = total_particles  # Set initial total count
        
        escape_fraction[timestep] = escaped_particles / total_particles_initial if total_particles_initial else 0
        
        if escape_time is None and escape_fraction[timestep] > 0.1:
            escape_time = timestep
    
    if escape_time is not None:
        print(f"Polymer escape time: {escape_time}")
    
    return escape_fraction

# Plot escape fraction over time
def save_escape_fraction_plot(escape_fraction, filename='escape_fraction.png'):
    times = sorted(escape_fraction.keys())
    fractions = [escape_fraction[t] for t in times]
    
    plt.figure(figsize=(8, 5))
    plt.plot(times, fractions, marker='o', linestyle='-', color='b')
    plt.xlabel('Time')
    plt.ylabel('Fraction of Escaped Particles')
    plt.title('Fraction of Escaped Particles Over Time')
    plt.grid()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    filename = 'bin_particles.lammpstrj'
    bin_division = 0.05
    xmin, xmax = 0, 53.8  # Define simulation box
    dx = bin_division * (xmax - xmin)
    
    time_data = parse_lammps_chunk_data(filename)

    escape_fraction = count_escaped_particles(time_data, dx)
    
    save_escape_fraction_plot(escape_fraction)

