import numpy as np

# Define your region's limits (adjust as needed)
x_min, x_max = 12, 18
y_min, y_max = 0.1, 9.9
z_min, z_max = 0.1, 9.9

total_rho = 3.0

# Function to count type 1 particles in the given region
def count_type_1_in_region(file_name, x_min, x_max, y_min, y_max, z_min, z_max):
    count = 0
    total_count = 0.0
    weight = 0.0
    with open(file_name, 'r') as file:
        # Read the dump file
        timestep = None
        i=0
        
        
        for line in file:
            
            if line.startswith("ITEM: TIMESTEP"):
                i = 0
                total_count = total_count + count
                weight  =  weight + 1
                count = 0
            if (i>9) :
                # Reading atom data
                data = line.strip().split()
                atom_type = int(data[1])
                
                x = float(data[2])
                
                # Check if the atom is type 1 and inside the defined region
                if atom_type == 1 and x_min <= x <= x_max :
                    count += 1
            i =i+1

    return total_count/weight

# Run the function with the file name and region ranges
file_name = "simulation_data.lammpstrj"  # Adjust to your actual file path
num_type_1_in_region = count_type_1_in_region(file_name, x_min, x_max, y_min, y_max, z_min, z_max)

volume = (x_max - x_min) * (y_max - y_min) * (z_max - z_min)

density =  num_type_1_in_region/volume

print("\n\nDensity in the selected region is given as :", density)

print("\n\nx is given as :", density/total_rho, "\n\n")
