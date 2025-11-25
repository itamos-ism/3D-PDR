# py3DPDR.py - python wrapper for 3D-PDR

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rc('font',size=17)

def makefile(
        F90="gfortran",
        CC="gcc",
        CPPFLAGS="-cpp",
        SUNDIALS=7.0,
        OPENMP=1,
        OPTIMISE=3,
        PYWRAP=1,
        DIMENSIONS=1,
        RAYTHEIA=0,
        RAYTHEIA_MO=0,
        XYZ=0,
        NETWORK="REDUCED",
        XRAYS=0,
        DUST="HTT91",
        GUESS_TEMP=1,
        THERMALBALANCE=1,
        FORCECONVERGENCE=1,
        GRAINRECOMB=0,
        SUPRATHERMAL=0,
        H2FORM="CT02",
        CRATTENUATION=0,
        RESTART=0
        ):

    config_values = {
        "F90": F90,
        "CC": CC,
        "CPPFLAGS": CPPFLAGS,
        "SUNDIALS": SUNDIALS,
        "OPENMP": OPENMP,
        "OPTIMISE": OPTIMISE,
        "PYWRAP": PYWRAP,
        "DIMENSIONS": DIMENSIONS,
        "RAYTHEIA": RAYTHEIA,
        "RAYTHEIA_MO": RAYTHEIA_MO,
        "XYZ": XYZ,
        "NETWORK": NETWORK,
        "XRAYS": XRAYS,
        "DUST": DUST,
        "GUESS_TEMP": GUESS_TEMP,
        "THERMALBALANCE": THERMALBALANCE,
        "FORCECONVERGENCE": FORCECONVERGENCE,
        "GRAINRECOMB": GRAINRECOMB,
        "SUPRATHERMAL": SUPRATHERMAL,
        "H2FORM": H2FORM,
        "CRATTENUATION": CRATTENUATION,
        "RESTART": RESTART
    }

    print(f"SUNDIALS = {SUNDIALS}")
    print(f"NETWORK = {NETWORK}")
    print(f"DUST = {DUST}")
    print(f"THERMALBALANCE = {THERMALBALANCE}")
    print(f"GRAINRECOMB = {GRAINRECOMB}")
    print(f"SUPRATHERMAL = {SUPRATHERMAL}")
    print(f"H2FORM = {H2FORM}")
    print(f"CRATTENUATION = {CRATTENUATION}")

    directory = "src/"
    # Write the key-value pairs to config.mk
    with open(directory + "config.mk", "w") as f:
        f.write(f"#==================================================================\n")
        f.write(f"#   	 	       3D-PDR Makefile\n")
        f.write(f"#----------------Compiler options----------------------------------\n")
        f.write(f"#SUNDIALS 	  : 7.0 new version\n")
        f.write(f"#		    2.5 used in Bisbas+12\n")
        f.write(f"#OPENMP   	  : 1 - run in parallel (OpenMP)\n")
        f.write(f"#         	    0 - run in serial (1 CPU) \n")
        f.write(f"#OPTIMISE 	  : use 0 or 1 for development & debugging\n")
        f.write(f"#     	            use 2 or 3 for running models\n")
        f.write(f"#PYWRAP		  : 1 - output for python wrapper\n")
        f.write(f"#		  : 0 - switch off\n")
        f.write(f"#-----------------Chemistry options--------------------------------\n")
        f.write(f"#DIMENSIONS 	  : 1 - compile for 1D models (needs RAYTHEIA = 0)\n")
        f.write(f"#		    3 - compile for 3D models\n")
        f.write(f"#RAYTHEIA 	  : 1 - switch on\n")
        f.write(f"#	 	  : 0 - switch off (uses the 2012 method)\n")
        f.write(f"#RAYTHEIA_MO      : 1 - memory optimization (code is slower, less RAM)\n")
        f.write(f"#                   0 - switch off (code is faster, more RAM)\n")
        f.write(f"#XYZ              : 0 - z coordinate changes first - xyz\n")
        f.write(f"#	   	  : 1 - x coordinate changes first - xyz\n")
        f.write(f"#NETWORK 	  : REDUCED (33 species & 331 reactions)\n")
        f.write(f"#         	    MEDIUM  (77 species & 1158 reactions; due to BALG)\n")
        f.write(f"#         	    FULL    (215 species & 2926 rections)\n")
        f.write(f"#XRAYS  	  : 1 - not yet fully working and works only with REDUCED\n")
        f.write(f"#		    0 - switch off\n")
        f.write(f"#DUST   	  : HTT91 - Hollenbach, Takahashi & Tielens (1991)\n")
        f.write(f"#                   0     - for isothermal dust models\n")
        f.write(f"#THERMALBALANCE   : 1 - for thermal balance iterations\n")
        f.write(f"#		    0 - for isothermal runs (needs GUESS_TEMP=0 also)\n")
        f.write(f"#FORCECONVERGENCE : 1 - Helps fast convergence (recommended)\n")
        f.write(f"#		    0 - Switch off\n")
        f.write(f"#GRAINRECOMB      : 1 - Electron recombination on dust grains\n")
        f.write(f"#	            0 - Switch off\n")
        f.write(f"#SUPRATHERMAL     : 1 - Suprathermal formation of CO via CH+\n")
        f.write(f"#	            0 - Switch off\n")
        f.write(f"#H2FORM	          : CT02   - Cazaux & Tielens (2002,2004) treatment\n")
        f.write(f"#	            SIMPLE - 3e-18*sqrt(Tgas)*exp(-Tgas/1e3) expression\n")
        f.write(f"#	            R07    - Roellig+07 benchmarking paper\n")
        f.write(f"#CRATTENUATION    : 1 - models CR attenuation (L/H of Padovani+ models)\n")
        f.write("#                   0 - switch off (default)\n")
        f.write(f"#RESTART          : 1 - restarts an interrupted model\n")
        f.write("#                   0 - switch off (default)\n")
        f.write(f"#------------------------------------------------------------------\n")

        for key, value in config_values.items():
            f.write(f"{key:<20} = {value}\n")  # Align formatting for readability

    print(f"Compiling 3D-PDR")
    currentdirectory = os.getcwd()
    os.chdir('src/')
    os.system('make clean;make')
    os.chdir(currentdirectory)
    print(f"Code compiled")

def params(
    # Input/Output - Densities
    icsdir="ics",
    icsfile="1Dn30.dat",
    outdir="sims",
    prefix="model",
    # PDR parameters
    fuv=1,
    zcr=1.0E-17,
    xrays=0.0,
    d2g=1.0,
    vturb=1.0,
    tend=1e7,
    grainrad=1.0E-5,
    avfac=6.289e-22,
    uvfac=3.02,
    redshift=0,
    av_crit=0.7,
    v_alfv=3.3,
    # Ray-tracing parameters
    healpix_level=0,
    theta_crit=1.3,
    # ODE Solver parameters
    rel_tol=1.0E-8,
    abs_tol=1.0E-30,
    # Chemistry iterations
    init_iter=8,
    max_iter=6000,
    # Thermal balance values
    tgas=40.0,
    tfloor=2.725,
    tgas_max=30000.0,
    tdust=20.0,
    Fcrit=0.005,
    Tdiff=0.01,
    # Coolant files
    coolant_files=None
    ):

    if coolant_files is None:
        coolant_files = []

    with open("params.dat", "w") as file:
        # Input/Output - Densities
        file.write("=================================\n")
        file.write("Input/Output - Densities\n")
        file.write("=================================\n")
        file.write(f"{icsdir}\n")
        file.write(f"{icsfile}\n")
        file.write(f"{outdir}\n")
        file.write(f"{prefix}\n")

        # PDR parameters
        file.write("=================================\n")
        file.write("PDR parameters\n")
        file.write("=================================\n")
        file.write(f"{fuv:.1f}\n")
        file.write(f"{zcr:.2E}\n")
        file.write(f"{xrays:.1f}\n")
        file.write(f"{d2g:.1f}\n")
        file.write(f"{vturb:.1f}\n")
        file.write(f"{tend:.1E}\n")
        file.write(f"{grainrad:.1E}\n")
        file.write(f"{avfac:.2E}\n")
        file.write(f"{uvfac:.2f}\n")
        file.write(f"{redshift:.1f}\n")
        file.write(f"{av_crit:.1f}\n")
        file.write(f"{v_alfv:.1f}\n")

        # Ray-tracing parameters
        file.write("=================================\n")
        file.write("Ray-tracing parameters\n")
        file.write("=================================\n")
        file.write(f"{healpix_level}\n")
        file.write(f"{theta_crit:.1f}\n")

        # ODE Solver parameters
        file.write("=================================\n")
        file.write("ODE Solver parameters\n")
        file.write("=================================\n")
        file.write(f"{rel_tol:.1E}\n")
        file.write(f"{abs_tol:.1E}\n")

        # Chemistry iterations
        file.write("=================================\n")
        file.write("Chemistry iterations\n")
        file.write("=================================\n")
        file.write(f"{init_iter}\n")
        file.write(f"{max_iter}\n")

        # Thermal balance values
        file.write("=================================\n")
        file.write("Thermal balance values\n")
        file.write("=================================\n")
        file.write(f"{tgas:.1f}\n")
        file.write(f"{tfloor:.1f}\n")
        file.write(f"{tgas_max:.1f}\n")
        file.write(f"{tdust:.1f}\n")
        file.write(f"{Fcrit:.3f}\n")
        file.write(f"{Tdiff:.2f}\n")

        # Coolant files
        file.write("=================================\n")
        file.write("Coolant files\n")
        file.write("=================================\n")
        file.write("12co.dat\n")
        file.write("12c+.dat\n")
        file.write("12c.dat\n")
        file.write("16o.dat\n")
        for coolant_file in coolant_files:
            file.write(f"{coolant_file}\n")
    
    print(f"File [params.dat] updated")

def run3DPDR():
    print(f"Running 3D-PDR. Check terminal for status")
    os.system('./3DPDR')
    print(f"Simulation finished")

def ics(
        nH=1000, 
        Avmax=10, 
        Avmin=1e-3, 
        resolution=30, 
        avfac=6.289e-22, 
        directory="ics/",
        filename="outgrid.dat"
        ):
    # Constants
    cm = 3.0856e18  # cm in a parsec
    
    # Compute the length in parsecs
    length = Avmax / (cm * avfac * nH)

    # Compute the number of logarithmic steps
    lambda_max = (np.log10(Avmax) - np.log10(Avmin)) * resolution
    lambda_max = int(lambda_max)

    # Open file and write initial values
    with open(directory+filename, "w") as f:
        f.write(f"{0.0:.7E} {0.0:.7E} {0.0:.7E} {nH:.7E}\n")
        
        for i in range(lambda_max + 1):
            av = 10 ** (np.log10(Avmin) + i / resolution)
            x = av / (cm * avfac * nH)
            f.write(f"{x:.7E} {0.0:.7E} {0.0:.7E} {nH:.7E}\n")

    print(f"File [{directory}{filename}] created")

def species(
        network=None, 
        Mg_plus=2.70E-07, 
        C_plus=1.40E-04, 
        He=1.00E-01, 
        O=3.00E-04, 
        H2=3.00E-01, 
        H=4.00E-01, 
        N_plus=6.76E-05,
        Na_plus=1.73E-06,
        F_plus=3.63E-08, 
        P_plus=2.57E-07, 
        S_plus=1.32E-05, 
        Fe_plus=3.16E-05, 
        Si_plus=3.23E-05, 
        Cl_plus=3.16E-07
        ):

    # Ensure Network is provided
    if network is None:
        raise ValueError("Error: 'network' must be specified (allowed values: REDUCED, MEDIUM, FULL)')")

    if network == "REDUCED":
        filename="chemfiles/species_reduced.d"
    if network == "MEDIUM":
        filename="chemfiles/species_medium.d"
    if network == "FULL":
        filename="chemfiles/species_full.d"
    while network not in ["REDUCED", "MEDIUM", "FULL"]:
        raise ValueError("Unknown 'network' (allowed values: REDUCED, MEDIUM, FULL)")

    # Species to modify
    species_to_modify = {
        "Mg+": Mg_plus,
        "C+": C_plus,
        "He": He,
        "O": O,
        "H2": H2,
        "H": H,
        "N": N_plus,
        "Na": Na_plus,
        "F+": F_plus,
        "P+": P_plus,
        "S+": S_plus,
        "Fe+": Fe_plus,
        "Si+": Si_plus,
        "Cl+": Cl_plus
    }

    # Read the file content
    if not os.path.exists(filename):
        print(f"Error: File '{filename}' not found!")
        return

    with open(filename, "r") as f:
        lines = f.readlines()

    # Modify relevant species
    modified_lines = []
    for line in lines:
        parts = line.strip().split(",")  # Split the line by comma
        if len(parts) < 4:
            modified_lines.append(line)  # Skip malformed lines
            continue

        index, species, value, mass = parts
        if species in species_to_modify:
            new_value = f"{species_to_modify[species]:.2E}"  # Format as scientific notation
            modified_line = f"{index},{species},{new_value},{mass}\n"
        else:
            modified_line = line  # Keep the line unchanged

        modified_lines.append(modified_line)

    # Write the updated content back to the file
    with open(filename, "w") as f:
        f.writelines(modified_lines)

    print(f"File '{filename}' successfully updated!")


def plot(x_var, y_vars, linestyle='-', scale="loglog", directory="sims/", prefix=None):


    if prefix is None:
        raise ValueError("Error: Give model input to plot specified in 'prefix'")

    filename = directory + prefix + ".pywrap"
    
    try:
        with open(filename, "r") as f:
            _model_header = f.readline().strip().split()  # Read first line as header
    except FileNotFoundError:
        print(f"Warning: '{filename}' not found! Load a valid file before plotting.")
        _model_header = None

    # Ensure y_vars is a list
    if isinstance(y_vars, str):
        y_vars = [y_vars]

    # Check if requested variables exist in the header
    missing_vars = [var for var in [x_var] + y_vars if var not in _model_header]
    if missing_vars:
        print(f"Error: Variables {missing_vars} not found in {filename}")
        return

    # Read the data
    data = np.loadtxt(filename, skiprows=1)  # Skip the header line
    col_indices = {name: i for i, name in enumerate(_model_header)}

    # Extract x and y values
    x_data = data[:, col_indices[x_var]]
    y_data = [data[:, col_indices[var]] for var in y_vars]

    # Plot the data with the selected scale
    plt.figure(figsize=(8, 6))
    
    for i, var in enumerate(y_vars):
        if scale == "loglog":
            plt.loglog(x_data, y_data[i], linestyle, label=var, linewidth=2)
        elif scale == "semilogy":
            plt.semilogy(x_data, y_data[i], linestyle, label=var, linewidth=2)
        elif scale == "semilogx":
            plt.semilogx(x_data, y_data[i], linestyle, label=var, linewidth=2)
        else:  # Default to linear if an invalid scale is provided
            plt.plot(x_data, y_data[i], linestyle, label=var, linewidth=2)

    # Labels and legend
    plt.xlabel(x_var)
    plt.ylabel(", ".join(y_vars))
    plt.legend()
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.title(f"Model = {filename}")
    plt.show()

