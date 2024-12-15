import csv

# Input and output file names
input_file = "msd.dat"
output_file = "msd_processed.txt"



time_step_m=0.005
# Open the input and output files
with open(input_file, "r") as infile, open(output_file, "w", newline="") as outfile:
    # Write header for the output CSV
    
    lines = infile.readlines()
    timestep = None
    msd_x, msd_y, msd_z, msd_total = None, None, None, None

    for line in lines:
        # Skip comments
        if line.startswith("#"):
            continue
        
        
        
        row_index, msd_value = map(float, line.split())
        
        if (msd_value==4.0):
            timestep = row_index 
        else :
        
            # Assign MSD components based on the row index
            if int(row_index) == 1:
                msd_x = msd_value
            elif int(row_index) == 2:
                msd_y = msd_value
            elif int(row_index) == 3:
                msd_z = msd_value
            elif int(row_index) == 4:
                msd_total = msd_value
                # Write the complete row for this timestep
                outfile.write(f"{timestep*time_step_m} {msd_x} {msd_y} {msd_z} {msd_total}\n")

