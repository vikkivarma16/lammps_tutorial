import numpy as np

def generate_hexagonal_lattice(lattice_constant, center, length, width, opening_radius):
    """
    Generates a 2D hexagonal lattice with normal along the x-axis and an opening at the center,
    formatted as a molecular template for LAMMPS.
    
    Parameters:
        lattice_constant (float): Distance between nearest neighbors.
        center (tuple): (x, y, z) coordinates of the center of mass.
        length (float): Extension along the y-axis.
        width (float): Extension along the z-axis.
        opening_radius (float): Radius of the circular opening at the center.
    
    Returns:
        list: List of (id, x, y, z, type) for lattice points.
    """
    a = lattice_constant
    lattice_points_hole = []
    lattice_points_solid = []
    atom_id = 1
    atom_id_solid  = 1 
    atom_type = 1  # Default atom type
    
    # Compute number of rows and columns needed
    num_rows = int(length / (np.cos(np.pi/6)* a))
    num_cols = int(width / a)
    
    x_offset, y_offset, z_offset = center
    
    for row in range(num_rows):
        for col in range(num_cols):
            y = col * a - width / 2 + y_offset
            z = row * np.cos(np.pi/6) * a - length / 2 + z_offset
            
            # Shift every alternate row
            if row % 2 == 1:
                y += a / 2
            
            # Check if the point is inside the opening
            if np.sqrt(y**2 + z**2) >= opening_radius:
                lattice_points_hole.append((atom_id, x_offset, y, z, atom_type))
                atom_id += 1
            lattice_points_solid.append((atom_id_solid, x_offset, y, z, atom_type))
            atom_id_solid = atom_id_solid + 1
    
    return lattice_points_hole, lattice_points_solid

# Example usage
lattice_constant = 0.4  # Set the lattice constant
yz_center = (0, 0, 0)  # Center of mass at (x, y, z)
length = 30  # Length along y-direction
width = 30   # Width along z-direction
opening_radius = 2.0  # Radius of the opening at the center

points, solids = generate_hexagonal_lattice(lattice_constant, yz_center, length, width, opening_radius)

# Save to a molecular template file for LAMMPS
with open("hexagonal_lattice_hollow.mol", "w") as f:
    f.write(f"\n")
    f.write(f"{3*len(points)} atoms\n")
    f.write(f"{2*len(points)} bonds\n\n")
    f.write("Coords\n\n")
    for p in points:
        f.write(f"{p[0]} {p[1]:.4f} {p[2]:.4f} {p[3]:.4f}\n")
        
    for p in points:
        f.write(f"{p[0]+len(points)} {(p[1]+0.5):.4f} {p[2]:.4f} {p[3]:.4f}\n")
        
    for p in points:
        f.write(f"{p[0]+2*len(points)} {(p[1]-0.5):.4f} {p[2]:.4f} {p[3]:.4f}\n")
    
    
    f.write("\nTypes\n\n")
    for p in points:
        f.write(f"{p[0]} {p[4]}\n")
        
    for p in points:
        f.write(f"{p[0]+len(points)} {p[4]+1}\n")
    
    for p in points:
        f.write(f"{p[0]+2*len(points)} {p[4]+1}\n")
    
    f.write("\n")
    
    f.write("Bonds \n\n")
    
    for p in points:
        f.write(f"{p[0]} {p[4]} {p[0]}  {p[0]+len(points)} \n")
        
    for p in points:
        f.write(f"{len(points)+p[0]} {p[4]} {p[0]}  {p[0]+2*len(points)} \n")
    
    
        
with open("hexagonal_lattice_hollow.txt", "w") as f:
    for p in points:
        f.write(f" {p[1]:.4f} {p[2]:.4f} {p[3]:.4f}\n")
    
    f.write("\nTypes\n\n")       
        




with open("hexagonal_lattice_solid.mol", "w") as f:
    f.write(f"\n")
    f.write(f"{3*len(solids)} atoms\n")
    f.write(f"{2*len(solids)} bonds\n\n")

    f.write("Coords\n\n")
    for p in solids:
        f.write(f"{p[0]} {p[1]:.4f} {p[2]:.4f} {p[3]:.4f}\n")
        
        
    for p in solids:
        f.write(f"{p[0]+len(solids)} {(p[1]+0.5):.4f} {p[2]:.4f} {p[3]:.4f}\n")
        
    for p in solids:
        f.write(f"{p[0]+2*len(solids)} {(p[1]-0.5):.4f} {p[2]:.4f} {p[3]:.4f}\n")
    
    
    f.write("\nTypes\n\n")
    for p in solids:
        f.write(f"{p[0]} {p[4]}\n")
    
    for p in solids:
        f.write(f"{p[0]+len(solids)} {p[4]+1}\n")
    
    for p in solids:
        f.write(f"{p[0]+2*len(solids)} {p[4]+1}\n")
        
    f.write("\n")
    
    f.write("Bonds \n\n")
    
    for p in solids:
        f.write(f"{p[0]} {p[4]} {p[0]}  {p[0]+len(solids)} \n")
        
    for p in solids:
        f.write(f"{len(solids)+p[0]} {p[4]} {p[0]}  {p[0]+2*len(solids)} \n")
    
        
        
with open("hexagonal_lattice_solid.txt", "w") as f:
    for p in solids:
        f.write(f" {p[1]:.4f} {p[2]:.4f} {p[3]:.4f}\n")
    
    f.write("\nTypes\n\n")


print(f"Generated {len(points)} lattice points with an opening radius of {opening_radius} for LAMMPS.")

