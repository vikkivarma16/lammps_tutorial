# density functional minimizer/potential generator in r and k space

# this code generates potential values for the distance measured from the center of a particle for the visualization and nature of the potential purpose both in the r and k space within the cartesian framework....
# in K space its a repetitive values from all the particles so its values are for the set of frequencies from zero to maximum relevant values... 


def pair_potential_particles_visualization():

    EPSILON = 0.0001
    
    import numpy as np
    from scipy.fft import fftfreq
    import matplotlib.pyplot as plt
    import json
    from scipy.integrate import simpson
    from scipy.special import jv  # Bessel function
    import calculator_pair_potential_custom
    # Define interaction potentials for different interaction types

    # Lennard-Jones potential
    def lennard_jones_potential(r, epsilon, sigma):
        return 4 * epsilon * ((sigma / r) ** 12 - (sigma / r) ** 6)

    # Gaussian potential
    def gaussian_potential(r, epsilon, sigma):
        return epsilon * np.exp(-(r ** 2) / (2 * sigma ** 2))

    def hard_core_potential(r, epsilon, sigma):
        
        large_value=1e6
        
        smoothing_width=0.05
        
        potential = np.where(r < sigma - smoothing_width, large_value, 0.0)  # Flat regions
        # Smooth transition near r = sigma
        smoothing_region = (sigma - smoothing_width <= r) & (r < sigma)
        potential[smoothing_region] = large_value * np.exp(-(r[smoothing_region] - (sigma - smoothing_width)) / smoothing_width)
        return potential

    def wca_potential_1(r, epsilon, sigma):
        
        return 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6) + epsilon

    def wca_potential(r, epsilon, sigma):
        
        # Define the cutoff for the WCA potential
        r_cutoff = 2**(1/6) * sigma
        r_upper_bound = 5 * sigma

        # Initialize the potential array
        V_wca = np.zeros_like(r)

        # Calculate the WCA potential based on the distance
        V_wca[r < r_cutoff] = -epsilon  # Region 1: r < 2^(1/6) * sigma
        V_wca[(r >= r_cutoff) & (r < r_upper_bound)] = 4 * epsilon * ((sigma / r[(r >= r_cutoff) & (r < r_upper_bound)])**12 - 2 * (sigma / r[(r >= r_cutoff) & (r < r_upper_bound)])**6)  # Region 2: 2^(1/6) * sigma < r < 5*sigma
        V_wca[r >= r_upper_bound] = 0  # Region 3: r > 5 * sigma

        return V_wca


        
        
        
        
        
        
        
        
        
        
    def bessel_fourier_transform(r, V_r):
        k_space = np.linspace(0.1, 10, len(r))  # Choose an appropriate k-space range
        V_k = np.zeros_like(k_space)

        # Calculate the Hankel transform using numerical integration
        for i, k in enumerate(k_space):
            integrand = r * V_r * jv(0, k * r)  # Bessel function of the first kind
            #V_k[i] = simpson(integrand, r)  # Use Simpson's rule for integration
            V_k[i] = simpson(y=integrand, x=r)

        return V_k, k_space











    # Function to read particle interaction data from the input_potential_generator.json file
    def read_particle_data_from_json_primary(json_file, interaction_types, closest_distances, interaction_strength, cutoff_ranges):
        with open(json_file, 'r') as file:
            data = json.load(file)
        data = data["particles_interactions_parameters"]
        
        interaction_types["primary"] = {}
        closest_distances["primary"] = {}
        interaction_strength["primary"] = {}
        cutoff_ranges["primary"] = {}

        # Loop through the interactions and parse the data
        for pair_type, values in data["interactions"]["primary"].items():
            interaction_types["primary"][pair_type] = values["type"]
            closest_distances["primary"][pair_type] = values["sigma"]
            interaction_strength["primary"][pair_type] = values["epsilon"]
            cutoff_ranges["primary"][pair_type] = values["cutoff"]
        
        file.close()
        return interaction_types["primary"], closest_distances["primary"], interaction_strength["primary"], cutoff_ranges["primary"]
        
    def read_particle_data_from_json_secondary(json_file, interaction_types, closest_distances, interaction_strength, cutoff_ranges):
        with open(json_file, 'r') as file:
            data = json.load(file)
        data = data["particles_interactions_parameters"]
        
        interaction_types["secondary"] = {}
        closest_distances["secondary"] = {}
        interaction_strength["secondary"] = {}
        cutoff_ranges["secondary"] = {}

        # Loop through the interactions and parse the data
        for pair_type, values in data["interactions"]["secondary"].items():
            interaction_types["secondary"][pair_type] = values["type"]
            closest_distances["secondary"][pair_type] = values["sigma"]
            interaction_strength["secondary"][pair_type] = values["epsilon"]
            cutoff_ranges["secondary"][pair_type] = values["cutoff"]
        file.close()
        return interaction_types["secondary"], closest_distances["secondary"], interaction_strength["secondary"], cutoff_ranges["secondary"]
        
    def read_particle_data_from_json_tertiary(json_file, interaction_types, closest_distances, interaction_strength, cutoff_ranges):
        with open(json_file, 'r') as file:
            data = json.load(file)
        data = data["particles_interactions_parameters"]
        interaction_types["tertiary"] = {}
        closest_distances["tertiary"] = {}
        interaction_strength["tertiary"] = {}
        cutoff_ranges["tertiary"] = {}

        # Loop through the interactions and parse the data
        for pair_type, values in data["interactions"]["tertiary"].items():
            interaction_types["tertiary"][pair_type] = values["type"]
            closest_distances["tertiary"][pair_type] = values["sigma"]
            interaction_strength["tertiary"][pair_type] = values["epsilon"]
            cutoff_ranges["tertiary"][pair_type] = values["cutoff"]
        
        file.close()
        return interaction_types["tertiary"], closest_distances["tertiary"], interaction_strength["tertiary"], cutoff_ranges["tertiary"]
        
    








    # Function to calculate interaction potential based on type
    def calculate_interaction_potential(r, interaction_type, epsilon, sigma):
        if interaction_type == 'lj':
            return lennard_jones_potential(r, epsilon, sigma)
        elif interaction_type == 'gs':
            return gaussian_potential(r, epsilon, sigma)
        elif interaction_type == 'hc':
            return hard_core_potential(r, epsilon, sigma)
        elif interaction_type == "wca":
            return wca_potential(r, epsilon, sigma)
        elif "custom" in interaction_type:
            dummy_string="pair_potential_"+interaction_type
            
            try:
                func = getattr(calculator_pair_potential_custom, dummy_string)
                v_r = np.zeros_like(r)
                func(r, v_r, epsilon, sigma)
                return v_r
            except AttributeError:
                print(f"       Error:    Unknown interaction type: {interaction_type}")
                exit(0)
            return func(r, epsilon, sigma)
        else:
            print(f"       Error:    Unknown interaction type: {interaction_type}")
            exit(0)
           
           
           
    
    
    
    
    

    # Function to export potential data to a file
    def export_potential_to_file(file_name, r_space, V_r):
        with open(file_name, 'w') as f:
            for r, V in zip(r_space, V_r):
                f.write(f"{r:.5f} {V:.5f}\n")  # Save r and potential in row-stacked format

    # Main function to calculate the interaction in both r-space and k-space
    
    
    
    
    
    
    
    
    
    
    def calculate_interactions(json_file):
        
        # Store all potentials and distances for plotting
        all_r_space = {}
        all_r_space["primary"] = {}
        all_r_space["secondary"] = {}
        all_r_space["tertiary"] = {}
        
        all_V_r = {}
        all_V_r["primary"] = {}
        all_V_r["secondary"] = {}
        all_V_r["tertiary"] = {}
        
        all_k_space = {}
        all_k_space["primary"] = {}
        all_k_space["secondary"] = {}
        all_k_space["tertiary"] = {}
        
        
        all_V_k = {}
        all_V_k["primary"] = {}
        all_V_k["secondary"] = {}
        all_V_k["tertiary"] = {}
        
        interaction_types = {}
        closest_distances = {}
        interaction_strength = {} 
        cutoff_ranges = {}
       
       
       
        
        
        
        # Get particle interaction data from the input_potential_generator.json file for the primary interaction ...
        interaction_types["primary"], closest_distances["primary"], interaction_strength["primary"], cutoff_ranges["primary"] = read_particle_data_from_json_primary(json_file, interaction_types, closest_distances, interaction_strength, cutoff_ranges)
        
        
        # Loop over interaction types and calculate interaction potentials
        for pair_type in interaction_types["primary"].keys():
            # Get the interaction parameters
            
            
            interaction_model = interaction_types["primary"][pair_type]  # Interaction type (e.g., Lennard-Jones or Gaussian)
            sigma = closest_distances["primary"][pair_type]  # Closest distance for the interaction type
            epsilon = interaction_strength["primary"][pair_type]  # Interaction strength (epsilon)
            cutoff = cutoff_ranges["primary"][pair_type]  # Cutoff range
            
            
            if interaction_model == 'gs':
                r_start = 0.0
            elif interaction_model == 'lj':
                r_start = sigma*2.0**(1.0/6.0)
            elif interaction_model == 'hc':
                r_start = sigma - 0.5
            elif interaction_model == 'wca':
                r_start = 0.0
            elif "custom" in interaction_model:
                r_start = 0.0
            else:
                print( "Error: the potential has not been implemented in the system, please choose an appropriate one.")
                exit(0)

            r_space = np.linspace(r_start, cutoff, 1000)

            
        

            
            # Calculate interaction potential in r-space
            V_r = calculate_interaction_potential(r_space, interaction_model, epsilon, sigma)
            
            
            # Calculate the interaction in k-space using the Bessel Fourier transform
            V_k, k_space = bessel_fourier_transform(r_space, V_r)
            
            # Save for plotting later
            all_r_space["primary"][pair_type] = r_space
            all_V_r["primary"][pair_type] = V_r
            all_k_space["primary"][pair_type] = k_space
            all_V_k["primary"][pair_type] = V_k
            
            
            '''
            # Export the r-space potential to a file
            r_file_name = f"potential_{pair_type}_r.txt"
            export_potential_to_file(r_file_name, r_space, V_r)
            
            # Export the k-space potential to a file
            k_file_name = f"potential_{pair_type}_k.txt"
            export_potential_to_file(k_file_name, k_space, V_k)
            '''
        
        
        
        
        
        
         # the same procedure is being repeated for the secondary interaction ...
        interaction_types["secondary"], closest_distances["secondary"], interaction_strength["secondary"], cutoff_ranges["secondary"] = read_particle_data_from_json_secondary(json_file, interaction_types, closest_distances, interaction_strength, cutoff_ranges)
        
        
        # Loop over interaction types and calculate interaction potentials
        for pair_type in interaction_types["secondary"].keys():
            # Get the interaction parameters
            
            
            interaction_model = interaction_types["secondary"][pair_type]  # Interaction type (e.g., Lennard-Jones or Gaussian)
            sigma = closest_distances["secondary"][pair_type]  # Closest distance for the interaction type
            epsilon = interaction_strength["secondary"][pair_type]  # Interaction strength (epsilon)
            cutoff = cutoff_ranges["secondary"][pair_type]  # Cutoff range
            
            
            if interaction_model == 'gs':
                r_start = 0.0
            elif interaction_model == 'lj':
                r_start = sigma*2.0**(1.0/6.0)
            elif interaction_model == 'hc':
                r_start = sigma - 0.5
            elif interaction_model == 'wca':
                r_start = 0.0
            elif "custom" in interaction_model:
                r_start = 0.0
            else:
                print( "Error: the potential has not been implemented in the system, please choose an appropriate one.")
                exit(0)

            r_space = np.linspace(r_start, cutoff, 1000)

            
        

            
            # Calculate interaction potential in r-space
            V_r = calculate_interaction_potential(r_space, interaction_model, epsilon, sigma)
            
            # Calculate the interaction in k-space using the Bessel Fourier transform
            V_k, k_space = bessel_fourier_transform(r_space, V_r)
            
            # Save for plotting later
            all_r_space["secondary"][pair_type] = r_space
            all_V_r["secondary"][pair_type] = V_r
            all_k_space["secondary"][pair_type] = k_space
            all_V_k["secondary"][pair_type] = V_k
            
            
            '''
            # Export the r-space potential to a file
            r_file_name = f"potential_{pair_type}_r.txt"
            export_potential_to_file(r_file_name, r_space, V_r)
            
            # Export the k-space potential to a file
            k_file_name = f"potential_{pair_type}_k.txt"
            export_potential_to_file(k_file_name, k_space, V_k)
            '''
        
        
        
        
        
        
        
        # the same procedure is being repeated for the tertiary interaction ...
        interaction_types["tertiary"], closest_distances["tertiary"], interaction_strength["tertiary"], cutoff_ranges["tertiary"] = read_particle_data_from_json_tertiary(json_file, interaction_types, closest_distances, interaction_strength, cutoff_ranges)
        
        
        # Loop over interaction types and calculate interaction potentials
        for pair_type in interaction_types["tertiary"].keys():
            # Get the interaction parameters
            
            
            interaction_model = interaction_types["tertiary"][pair_type]  # Interaction type (e.g., Lennard-Jones or Gaussian)
            sigma = closest_distances["tertiary"][pair_type]  # Closest distance for the interaction type
            epsilon = interaction_strength["tertiary"][pair_type]  # Interaction strength (epsilon)
            cutoff = cutoff_ranges["tertiary"][pair_type]  # Cutoff range
            
            
            if interaction_model == 'gs':
                r_start = 0.0
            elif interaction_model == 'lj':
                r_start = sigma*2.0**(1.0/6.0)
            elif interaction_model == 'hc':
                r_start = sigma - 0.5
            elif interaction_model == 'wca':
                r_start = 0.0
            elif "custom" in interaction_model:
                r_start = 0.0
            else:
                print( "Error: the potential has not been implemented in the system, please choose an appropriate one.")
                exit(0)

            r_space = np.linspace(r_start, cutoff, 1000)

            
        

            
            # Calculate interaction potential in r-space
            V_r = calculate_interaction_potential(r_space, interaction_model, epsilon, sigma)
            
            # Calculate the interaction in k-space using the Bessel Fourier transform
            V_k, k_space = bessel_fourier_transform(r_space, V_r)
            
            # Save for plotting later
            all_r_space["tertiary"][pair_type] = r_space
            all_V_r["tertiary"][pair_type] = V_r
            all_k_space["tertiary"][pair_type] = k_space
            all_V_k["tertiary"][pair_type] = V_k
            
            
            '''
            # Export the r-space potential to a file
            r_file_name = f"potential_{pair_type}_r.txt"
            export_potential_to_file(r_file_name, r_space, V_r)
            
            # Export the k-space potential to a file
            k_file_name = f"potential_{pair_type}_k.txt"
            export_potential_to_file(k_file_name, k_space, V_k)
            '''
        
       
       
       
       
       
       
       
       
       
       
       
       
       
       
        
        # Define a list of line styles and colors
        line_styles = ['-', '--', '-.', ':']  # Solid, dashed, dash-dot, dotted
        colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']  # Blue, green, red, cyan, magenta, yellow, black

        # Plot all potentials in r-space and k-space side by side
        fig, axes = plt.subplots(1, 2, figsize=(19, 8), dpi=100)  # Set figure size to 1900x1080 pixels

        # Plot r-space potentials for the primary interaction
        for i, pair_type in enumerate(all_r_space["primary"].keys()):
            axes[0].plot(all_r_space["primary"][pair_type], all_V_r["primary"][pair_type], 
                          label=f'{pair_type} {"primary"} ({interaction_types["primary"][pair_type].capitalize()})', 
                          linewidth=2.5, 
                          linestyle=line_styles[i % len(line_styles)], 
                          color=colors[i % len(colors)])  # Cycle through line styles and colors
                          
                          
                          
        # Plot r-space potentials for the secondary interaction
        for i, pair_type in enumerate(all_r_space["secondary"].keys()):
            axes[0].plot(all_r_space["secondary"][pair_type], all_V_r["secondary"][pair_type], 
                          label=f'{pair_type} {"secondary"} ({interaction_types["secondary"][pair_type].capitalize()})', 
                          linewidth=2.5, 
                          linestyle=line_styles[i % len(line_styles)], 
                          color=colors[i % len(colors)])  # Cycle through line styles and colors
                          
                          
        # Plot r-space potentials for the tertiary interaction
        for i, pair_type in enumerate(all_r_space["tertiary"].keys()):
            axes[0].plot(all_r_space["tertiary"][pair_type], all_V_r["tertiary"][pair_type], 
                          label=f'{pair_type} {"tertiary"} ({interaction_types["tertiary"][pair_type].capitalize()})', 
                          linewidth=2.5, 
                          linestyle=line_styles[i % len(line_styles)], 
                          color=colors[i % len(colors)])  # Cycle through line styles and colors
                          
                          
                          
                          
                          

        axes[0].set_title('Potential V(r) for all interaction pairs')
        axes[0].set_xlabel('r')
        axes[0].set_ylabel('V(r)')
        axes[0].legend()
        axes[0].grid(True)
        axes[0].set_ylim(-10, 10)  # Set y-range for r-space plot

        # Plot k-space potentials for the primary interactions
        for i, pair_type in enumerate(all_k_space["primary"].keys()):
            axes[1].plot(all_k_space["primary"][pair_type], all_V_k["primary"][pair_type], 
                          label=f'{pair_type} {"primary"} ({interaction_types["primary"][pair_type].capitalize()})', 
                          linewidth=2.5, 
                          linestyle=line_styles[i % len(line_styles)], 
                          color=colors[i % len(colors)])  # Cycle through line styles and colors
        
        # Plot k-space potentials for the secondary interactions
        for i, pair_type in enumerate(all_k_space["secondary"].keys()):
            axes[1].plot(all_k_space["secondary"][pair_type], all_V_k["secondary"][pair_type], 
                          label=f'{pair_type} {"secondary"} ({interaction_types["secondary"][pair_type].capitalize()})', 
                          linewidth=2.5, 
                          linestyle=line_styles[i % len(line_styles)], 
                          color=colors[i % len(colors)])  # Cycle through line styles and colors
                          
        
        # Plot k-space potentials for the tertiary interactions
        for i, pair_type in enumerate(all_k_space["tertiary"].keys()):
            axes[1].plot(all_k_space["tertiary"][pair_type], all_V_k["tertiary"][pair_type], 
                          label=f'{pair_type} {"tertiary"} ({interaction_types["tertiary"][pair_type].capitalize()})', 
                          linewidth=2.5, 
                          linestyle=line_styles[i % len(line_styles)], 
                          color=colors[i % len(colors)])  # Cycle through line styles and colors
                          
                          
                          

        axes[1].set_title('Potential V(k) for all interaction pairs (Bessel Fourier Transform)')
        axes[1].set_xlabel('k')
        axes[1].set_ylabel('V(k)')
        axes[1].legend()
        axes[1].grid(True)
        axes[1].set_ylim(-10, 10)  # Set y-range for r-space plot
        

        plt.tight_layout()

        # Save the plot as a PNG file with high resolution
        plt.savefig('vis_interaction_potentials.png', dpi=300)  # Increase dpi for high resolution

        
        
        
        

        print("\n\n\n... all the potentials have been implemented successfully in the proper way, which could be observed in the exported png file ...\n\n\n")

        
        
        
        
      

        # Show the plot
        #plt.show()
        
    # Call the function to execute, using the input_potential_generator.json file
    json_file_path = 'input_data_particles_interactions_parameters.json'  # Make sure this file exists and follows the correct format
    calculate_interactions(json_file_path)

    return 0
    
#pair_potential_particles_visualization()
