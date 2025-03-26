# density functional minimizer/executor

# this is the main code which read and write the data in an executable format and then also run the program for the calculation...






# this part read the input file and print the json file for the further processing of the input file in the selected format...







# Density Functional Minimizer/Executor
# Main code that reads, writes, and processes input data to execute the program

# Density Functional Minimizer/Executor
import numpy as np
import json
import math
import matplotlib.pyplot as plt
# Assuming you have pynufft installed
from pynufft import NUFFT
from scipy.fft import fft , ifft
import pyfftw.interfaces as fftw





from generator_input_data_particles_interactions_parameters import data_exporter_particles_interactions_parameters as interactions_parameters
from generator_input_data_space_confinement_parameters import data_exporter_space_confinement_parameters as confinement_parameters
from generator_input_data_simulation_thermodynamic_parameters import data_exporter_simulation_thermodynamic_parameters as thermodynamic_parameters




input_file_name = 'executor_input.in'

try:
    interactions_parameters(input_file_name)
except Exception as e:
    print("Error exporting input data:", e)
    exit(0)
    
try:
    confinement_parameters(input_file_name)
except Exception as e:
    print("Error exporting input data:", e)
    exit(0)

try:
    thermodynamic_parameters(input_file_name)
except Exception as e:
    print("Error exporting input data:", e)
    exit(0)






from generator_pair_potential_particles_visualization import pair_potential_particles_visualization as visualizer
from generator_k_and_r_space_box import r_k_space as rks

try:
    visualizer()
except Exception as e:
    print("Error visualizing particles:", e)
    exit(0)
    
try:
    rks()
except Exception as e:
    print("Error generating k and r space:", e)
    exit(0)



#from generator_wall_potential_values_visualization import wall_potential_values_visualization as wp_values
from generator_bulk_rho_mue_r_space import bulk_rho_mue_r_space as brm
from calculator_FMT_weights_1d_cartesian  import fmt_weights_1d as fm_weights






# File paths





try: 
    brm()
except Exception as e:
    print("Error while generating the bulk rho and mue value ... \n")
    exit(0)

try:
    fm_weights()
except Exception as e:
    print ("Error while calculating FMT weights in k space ... \n")
    exit(0)



from generator_wall_potential_values_visualization import wall_potential_values_visualization as wp_values

try:
    wp_values()
except Exception as e:
    print("Error generating wall potential:", e)
    exit(0)









# Load space properties
json_file_interaction = "input_data_particles_interactions_parameters.json"
json_file_thermodynamics = "input_data_simulation_thermodynamic_parameters.json"
json_file_confinement = "input_data_space_confinement_parameters.json"
try:
    
    with open(json_file_thermodynamics, "r") as file:
        thermodynamics = json.load(file)["simulation_thermodynamic_parameters"]
    temperature, rho, iteration_max = thermodynamics["temperature"], thermodynamics["rho"], thermodynamics["iteration_max"]
    print(f"\n\n... thermodynamic properties loaded: Temperature = {temperature}, Rho = {rho}, Max Iterations = {iteration_max} ...\n\n")
except FileNotFoundError:
    print("Thermodynamic properties file not found.")
    exit()

try:
    with open(json_file_interaction, 'r') as file:
        data_interactions = json.load(file)["particles_interactions_parameters"]["interactions"]
    interaction_types = {}
    closest_distances = {}
    interaction_strength = {}
    cutoff_ranges = {}
    
    interaction_types["primary"] = {k: v["type"] for k, v in data_interactions["primary"].items()}
    closest_distances["primary"] = {k: v["sigma"] for k, v in data_interactions["primary"].items()}
    interaction_strength["primary"] = {k: v["epsilon"] for k, v in data_interactions["primary"].items()}
    cutoff_ranges["primary"] = {k: v["cutoff"] for k, v in data_interactions["primary"].items()}
    
    interaction_types["secondary"] = {k: v["type"] for k, v in data_interactions["secondary"].items()}
    closest_distances["secondary"] = {k: v["sigma"] for k, v in data_interactions["secondary"].items()}
    interaction_strength["secondary"] = {k: v["epsilon"] for k, v in data_interactions["secondary"].items()}
    cutoff_ranges["secondary"] = {k: v["cutoff"] for k, v in data_interactions["secondary"].items()}
    
    interaction_types["tertiary"] = {k: v["type"] for k, v in data_interactions["tertiary"].items()}
    closest_distances["tertiary"] = {k: v["sigma"] for k, v in data_interactions["tertiary"].items()}
    interaction_strength["tertiary"] = {k: v["epsilon"] for k, v in data_interactions["tertiary"].items()}
    cutoff_ranges["tertiary"] = {k: v["cutoff"] for k, v in data_interactions["tertiary"].items()}
    
    with open(json_file_interaction, 'r') as file:
        species = {k: v["rho_frac"] for k, v in json.load(file)["particles_interactions_parameters"]["species"].items()}
    print(f"\n... loaded {len(species)} species for simulation ...\n\n")
    
except FileNotFoundError:
    print("Interaction properties file not found.")
    exit()
   
   
try: 
    with open(json_file_confinement, 'r') as file:
        data = json.load(file)

    # Accessing data from the JSON
    space_properties = data["space_confinement_parameters"]["space_properties"]
    box_properties = data["space_confinement_parameters"]["box_properties"]

    # Print the loaded data
    print("Space Properties:")
    print(f"  Dimension: {space_properties['dimension']}")
    print(f"  Confinement Type: {space_properties['confinement']}")

    print("\nBox Properties:")
    print(f"  Box Lengths: {box_properties['box_length']}")
    print(f"  Box Points: {box_properties['box_points']}")

    # Example: Access individual values if needed
    dimension = space_properties["dimension"]
    confinement_type = space_properties["confinement"]
    box_lengths = box_properties["box_length"]
    box_points = box_properties["box_points"]

except FileNotFoundError:
    print("Space and properties file could not be found.")
    exit()


        



print("\n") 
interaction_level = "primary"
# Display interactions for verification
for pair_type, interaction_type in interaction_types[interaction_level].items():
    print(f"... primary interaction Pair {pair_type}:  Type = {interaction_type}, Sigma = {closest_distances[interaction_level][pair_type]}, "
          f"Epsilon = {interaction_strength[interaction_level][pair_type]}, Cutoff = {cutoff_ranges[interaction_level][pair_type]}\n")
          
print("\n")          
interaction_level = "secondary"
# Display interactions for verification
for pair_type, interaction_type in interaction_types[interaction_level].items():
    print(f"... secondary interaction Pair {pair_type}:  Type = {interaction_type}, Sigma = {closest_distances[interaction_level][pair_type]}, "
          f"Epsilon = {interaction_strength[interaction_level][pair_type]}, Cutoff = {cutoff_ranges[interaction_level][pair_type]}\n")
          
print("\n")       
interaction_level = "tertiary"
# Display interactions for verification
for pair_type, interaction_type in interaction_types[interaction_level].items():
    print(f"... tertiary interaction Pair {pair_type}:  Type = {interaction_type}, Sigma = {closest_distances[interaction_level][pair_type]}, "
          f"Epsilon = {interaction_strength[interaction_level][pair_type]}, Cutoff = {cutoff_ranges[interaction_level][pair_type]}  \n")
          
print("\n\n\n")




import numpy as np

nx= int(box_points[0])

# Initialize empty lists for each column
x, y, z, rho_r, mue_r = [], [], [], [], []
# Read the text file for the bulk rho mue values ... 
pid=0
for particle, rho in species.items():
    
    rho_ind= []
    mue_ind= []
    with open("supplied_data_bulk_mue_rho_r_space.txt", "r") as file:
        for line in file:
            # Skip comment lines
            if line.startswith("#"):
                continue
            # Split the line into columns and convert them to floats
            columns = line.strip().split()
            if (pid == 0 ):
                x.append(float(columns[0]))
                y.append(float(columns[1]))
                z.append(float(columns[2]))
            
            li=[]
            chi=[]
            i=3+pid*2
            rho_ind.append(float(columns[i]))
            mue_ind.append (float(columns[i+1]))
    rho_r.append(rho_ind)
    mue_r.append(mue_ind)
    pid = pid+2
    # Convert lists to numpy arrays
x = np.array(x)
y = np.array(y)
z = np.array(z)
rho_r = np.array(rho_r)
mue_r = np.array(mue_r)





# Read the text file for the fmt weights for the k space data 



kx, ky, kz = [], [], []  # Define kx, ky, kz before the loop
fmt_weights = {}  # Define fmt_weights to hold all species weights

for key, rho in species.items():
    fmt_weights_ind = []  # Initialize a list for individual species weights
    
    # Open file for the current species
    with open(f"supplied_data_weight_FMT_k_space_{key}.txt", "r") as file:
        for line in file:
            # Skip comment lines
            if line.startswith("#"):
                continue

            # Split the line into columns and convert them to floats
            columns = line.strip().split()
           
            
            # Collect rho-related values for this species
            li = []

            li.append(complex(columns[3]))
            li.append(complex(columns[4]))
            li.append(complex(columns[5]))
            li.append(complex(columns[6]))
            li.append(complex(columns[7]))
            li.append(complex(columns[8]))
       
            fmt_weights_ind.append(li)
            
    
    fmt_weights_ind = np.array (fmt_weights_ind)
    fmt_weights[key] = fmt_weights_ind # Append the individual weights list to fmt_weights




k_space_file_path = 'supplied_data_k_space.txt'
k_space = np.loadtxt(k_space_file_path)
k_space = np.array(k_space)

kx = k_space[:,0]
ky = k_space[:,1]
kz = k_space[:,2]





v_ext={}
for key in species:
    with open(f"supplied_data_walls_potential_{key}_r_space.txt", "r") as file:
        v_ind=[]
        for line in file:
            # Skip comment lines
            if line.startswith("#"):
                continue
            # Split the line into columns and convert them to floats
            columns = line.strip().split()
            v_ind.append( float(columns[3]))
        v_ext[key] = np.array(v_ind)
        
print ("\n\n... supplied data has been imported successfully ...\n\n\n")








# this is the main regions for the calculation for the 1 d walls confinement DFT simulation ...

i = 0
j = 0

iteration = 0
rho_r_initial = np.array(rho_r)
rho_r_current = np.array(rho_r)

piee = np.pi







# for mean field treatment
epsilonij = [[0.0]*len(species) for _ in range(len(species))]
sigmaij = [[1.0]*len(species) for _ in range(len(species))]
interaction_type_ij = [["gs"]*len(species) for _ in range(len(species))]


# for hard core treatment
sigmai_p = [0.0]*len(species)
flag = [0]*len(species)


rhos = []

i = 0 
j = 0
grand_rosenfeld_flag = 0

for species_type, rho_value in species.items():

    rhos.append(species[species_type])
    # Initialize total chemical potential for the current species
    j = 0
    for other_species, rho_other in species.items():
        # Determine the interaction data key
        interaction_key1 = f"{species_type}{other_species}"
        interaction_key2 = f"{other_species}{species_type}"

        if interaction_key1 in interactions["primary"]:
            interaction_data = interactions["primary"][interaction_key1]
        elif interaction_key2 in interactions:
            interaction_data = interactions["primary"][interaction_key2]
        else:
            # Skip if no interaction data is found
            continue

        # Unpack interaction parameters
        epsilon = interaction_data["epsilon"]
        sigma_ij = interaction_data["sigma"]
        interaction_type = interaction_data["type"]
        cutoff = interaction_data["cutoff"]

        # Check if it's a hard-core (hc) interaction
        if interaction_type == "hc":
            grand_rosenfeld_flag = 1
            if species_type == other_species:
            # Approximate hard-core chemical potential contribution
                sigmai_p[i] = sigma_ij
                flag[i] = 1
            else:
                print("\n ...wrong parameters have been defined... \n")
                exit(0)
        elif interaction_type == "ghc":
            grand_rosenfeld_flag = 1
            if species_type == other_species:
            # Approximate hard-core chemical potential contribution
                sigmai_p[i] = sigma_ij
            else:
                print("\n ...wrong parameters have been defined... \n")
                exit(0)
                
        else :
            interaction_type_ij[i][j] = interaction_type
            epsilonij[i][j] = epsilon
            sigmaij[i][j] = sigma_ij
        j = j+1

    j = 0
    for other_species, rho_other in species.items():
        # Determine the interaction data key
        
        interaction_key1 = f"{species_type}{other_species}"
        interaction_key2 = f"{other_species}{species_type}"

        if interaction_key1 in interactions["secondary"]:
            interaction_data = interactions["secondary"][interaction_key1]
        elif interaction_key2 in interactions:
            interaction_data = interactions["secondary"][interaction_key2]
        else:
            # Skip if no interaction data is found
            continue

        # Unpack interaction parameters
        epsilon = interaction_data["epsilon"]
        sigma_ij = interaction_data["sigma"]
        interaction_type = interaction_data["type"]
        cutoff = interaction_data["cutoff"]

        
        # Check if it's a hard-core (hc) interaction
        if interaction_type == "hc":
            grand_rosenfeld_flag = 1
            if species_type == other_species:
            # Approximate hard-core chemical potential contribution
                sigmai_p[i] = sigma_ij
                flag[i] = 1
            else:
                print("\n ...wrong parameters have been defined... \n")
                exit(0)
        elif interaction_type == "ghc":
            grand_rosenfeld_flag = 1
            if species_type == other_species:
            # Approximate hard-core chemical potential contribution
                sigmai_p[i] = sigma_ij
            else:
                print("\n ...wrong parameters have been defined... \n")
                exit(0)
                
        else :
            interaction_type_ij[i][j] = interaction_type
            epsilonij[i][j] = epsilon
            sigmaij[i][j] = sigma_ij
    
        j = j+1

            
    j = 0
    for other_species, rho_other in species.items():
        # Determine the interaction data key
        
        interaction_key1 = f"{species_type}{other_species}"
        interaction_key2 = f"{other_species}{species_type}"

        if interaction_key1 in interactions["tertiary"]:
            interaction_data = interactions["tertiary"][interaction_key1]
        elif interaction_key2 in interactions:
            interaction_data = interactions["tertiary"][interaction_key2]
        else:
            # Skip if no interaction data is found
            continue

        # Check if it's a hard-core (hc) interaction
        if interaction_type == "hc":
            grand_rosenfeld_flag = 1
            if species_type == other_species:
            # Approximate hard-core chemical potential contribution
                sigmai_p[i] = sigma_ij
                flag[i] = 1
            else:
                print("\n ...wrong parameters have been defined... \n")
                exit(0)
        elif interaction_type == "ghc":
            grand_rosenfeld_flag = 1
            if species_type == other_species:
            # Approximate hard-core chemical potential contribution
                sigmai_p[i] = sigma_ij
            else:
                print("\n ...wrong parameters have been defined... \n")
                exit(0)
                
        else :
            interaction_type_ij[i][j] = interaction_type
            epsilonij[i][j] = epsilon
            sigmaij[i][j] = sigma_ij
        j = j+1


    #print("mue before ideal", mue_total/temperature)
    # Store the total chemical potential for the species
    i = i+1
    



if (grand_rosenfeld_flag == 1):

    def hard_core_approach(sigmai, rhos, flag):
        
        variables = [[sp.symbols(f"n{j}_{i} ") for j in range(6)] for i in range(len(sigmai))]
        etas = [ variables[i][3] for i in range(len(sigmai))]
        

        fac1 = 1 - sum(etas)
        fac2 = sum(etas[i] for i in range(len(sigmai)) if flag[i] == 1)
        
        
        phi0 = fac1 * sp.log(1 - fac2) + fac2
        
        
        diff_1 = [sp.diff(phi0, etas[i]) for i in range(len(sigmai))]
        diff_2 = [[sp.diff(phi0, etas[i], etas[j]) for j in range(len(sigmai))] for i in range(len(sigmai))]
        diff_3 = [[[sp.diff(phi0, etas[i], etas[j], etas[k]) for k in range(len(sigmai))] for j in range(len(sigmai))] for i in range(len(sigmai))]

      
        
        

        phi1 = sum(variables[i][0] * diff_1[i] for i in range(len(sigmai)))
        phi2 = sum((variables[i][1] * variables[j][2] - variables[i][4] * variables[j][5]) * diff_2[i][j] for i in range(len(sigmai)) for j in range(len(sigmai)))
        phi3 = sum((1/(8*np.pi)) * (variables[i][2]*variables[j][2]*variables[k][2]/3.0  - variables[i][2] *variables[j][5] * variables[k][5] + (3/2) * ( variables[i][5] *variables[j][5] * ((variables[i][2] - 4 * variables[i][3]/sigmai[i]) - variables[k][2]/2 ) - ((variables[i][2] - 4 * variables[i][3]/sigmai[i]) - variables[i][2]/2) * ((variables[i][2] - 4 * variables[i][3]/sigmai[i]) - variables[j][2]/2) * ((variables[i][2] - 4 * variables[i][3]/sigmai[i]) - variables[k][2]/2 ) + 2 * ( ((variables[i][2] - 4 * variables[i][3]/sigmai[i]) - variables[i][2]/2 )/2 *((variables[i][2] - 4 * variables[i][3]/sigmai[i]) - variables[j][2]/2 )/2 *((variables[i][2] - 4 * variables[i][3]/sigmai[i]) - variables[k][2]/2 ))    )  ) * diff_3[i][j][k] for i in range(len(sigmai)) for j in range(len(sigmai)) for k in range(len(sigmai)))
        
        total_phi = phi1+phi2+phi3
        
        
        functions = [[sp.diff(total_phi, variables[i][j]) for j in range(6)] for i in range(len(sigmai))]
        
        functions_func = [[sp.lambdify(variables, functions[i][j], 'numpy') for j in range(6)] for i in range(len(sigmai))]

        return functions








    pdphi = hard_core_approach(sigmai, rhos, flag)












threshold = 0.001
alpha = 0.02

print("\n\n...number of iteration is given as:", iteration_max, "\n\n")

while (iteration < iteration_max):
    
    rho_r_initial = rho_r_current 
    
    # energy filtration for the hard core fmt terms ...
    rho_alpha_r={}
    fmt_flag = 0
    pid = 0
    
    if (grand_rosenfeld_functional == 1):
        pid = 0
        for particle, rho_other in species.items():
            rho_k_ind = fft(rho_r_current[pid])
            omega_rho_k = np.zeros((6, nx), dtype=complex)
            li=[]
            for i in range(6):
                omega_rho_k[i,:] = fmt_weights[particle][:, i] * rho_k_ind 
                
                    #print(fmt_weights[particle1][:, i])
               
            rho_alpha_r_ind= np.zeros((6,nx))
            for i in range(6):
                rho_alpha_r_ind[i,:]= ifft(omega_rho_k[i, :]).real
    
                
            rho_alpha_r[particle] = np.array(rho_alpha_r_ind)
            pid = pid + 1
        
        tdphi = []
        dphi=np.zeros((6, nx))
        
        pid = 0
        for particle, rho_other in species.items():
            dphi = np.zeros((6, nx))
            tdphi.append(dphi)
    
        for k in range (nx):
            pid = 0
            for particle, rho_other in species.items():
                if ( (1.0 - rho_alpha_r[particle][3, k]) > 0.00000001 and (rho_alpha_r[particle][3, k]) > 0.00000001):
                    variable  = []
                    for particle_in, rho_other_in in species.items():
                        variable.append( [rho_alpha_r[particle][0, k], rho_alpha_r[particle][1, k], rho_alpha_r[particle][2, k], rho_alpha_r[particle][3, k], rho_alpha_r[particle][4, k], rho_alpha_r[particle][5, k]] )
                    tdphi[pid][0, k] =  pdphi[pid][0](*variable)
                    tdphi[pid][1, k] =  pdphi[pid][1](*variable)
                    tdphi[pid][2, k] =  pdphi[pid][2](*variable)
                    tdphi[pid][3, k] =  pdphi[pid][3](*variable)
                    tdphi[pid][4, k] =  pdphi[pid][4](*variable)
                    tdphi[pid][5, k] =  pdphi[pid][5](*variable)
                pid = pid + 1
        
        pid = 0
        dphi_k_new = [] 
        for particle, rho_other in species.items():
            dphi_k = []
            for i in range(6):
                dphi_k_alpha = fft(tdphi[pid][i])
                dphi_k.append(dphi_k_alpha)
            dphi_k = np.array(dphi_k)
            dphi_k_new.append(dphi_k)
                

    
        pid = 0
        
        f_ext_frag = []
        for particle, rho_other in species.items():
            dphi_k_ind = dphi_k_new[pid]
            omega_dphi_k = np.zeros((6, nx), dtype=complex)
            li=[]
            for i in range(6):
                omega_dphi_k[i,:] = fmt_weights[particle][:, i] * dphi_k_ind[i] 
                
                    #print(fmt_weights[particle1][:, i])
               
            dphi_rho_alpha_ind= np.zeros((6,nx))
            for i in range(6):
                dphi_rho_alpha_ind[i,:]= ifft(omega_dphi_k[i, :]).real
    
                
            frag_energy = np.array(dphi_rho_alpha_ind)
            f_ext_frag.append(frag_energy)
            pid = pid + 1
  
        
        pid = 0
        
        total_df_ext = []
        for particle, rho_other in species.items():
            
            df_ext_ind = np.zeros(nx)
            
            for i in range(nx):
                df_ext_ind[i] = np.sum(fext_frag[pid][:, i])
            
            total_df_ext.append(df_ext_ind)
            pid = pid +1
        
      
        print(total_df_ext[0][0])
        
        exit(0)
    
    


print(rho_r_current)

line_styles = ['-', '--', '-.', ':']
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

# Plotting
plt.figure(figsize=(12, 8), dpi=300)  # High-definition plot

# Loop through each potential profile
for i, key in enumerate(species):
    # Cycle through line styles and colors
    style = line_styles[i % len(line_styles)]
    color = colors[i % len(colors)]
    plt.plot(x, rho_r_current[i], marker='o', linestyle=style, color=color, label=f'Species {key}')
    #plt.ylim(0, 0.004)

# Plot customization
plt.xlabel('Position Magnitude')
plt.ylabel('Density distribution')
plt.title('Density distribution for different species')
plt.grid(True)
plt.legend()

# Save the plot in high resolution
plt.savefig('vis_rho_distribution.png')












    
    
# testing region ... 
'''
import numpy as np
import matplotlib.pyplot as plt

# Example r_space values from 0 to a certain cutoff (e.g., 10)
r_space = np.linspace(0, 10, 1000)  # 1000 points in r-space
r_space_data = np.exp(-r_space)  # Sample data, e.g., a decaying exponential

# Perform FFT
k_space_data = np.fft.fft(r_space_data)

# Get the k-space values, assuming spacing between r-space points is uniform
dr = r_space[1] - r_space[0]  # calculate spacing in r-space
k_values = np.fft.fftfreq(len(r_space), d=dr) * 2 * np.pi  # Convert to angular frequency



plt.close('all')
plt.clf()  # Clears current figure for a new plot
# Plot r-space data
plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
plt.plot(r_space, r_space_data, label='r-space data')
plt.xlabel('r')
plt.ylabel('Amplitude')
plt.title('Real-space Data')
plt.legend()

# Plot k-space data (magnitude)
plt.subplot(1, 2, 2)
plt.plot(np.fft.fftshift(k_values), np.fft.fftshift(np.abs(k_space_data)), label='k-space data')
plt.xlabel('k')
plt.ylabel('Magnitude')
plt.title('k-space Data')
plt.legend()

plt.tight_layout()
plt.show()
'''
