import numpy as np
import h5py

def calculate_column_density(density_grid, abundance_grid=None, axis=0, cell_size=1.0):
    """
    Calculate the column density along a specified axis in a uniform 3D grid.
    Optionally multiply by abundance before summing.
    
    Parameters:
    -----------
    density_grid : numpy.ndarray
        3D array containing the density values in each cell
    abundance_grid : numpy.ndarray or None
        3D array containing abundance values (same shape as density_grid)
        If None, only density is used
    axis : int (0, 1, or 2)
        Axis along which to calculate the column density
    cell_size : float
        Physical size of each grid cell
    
    Returns:
    --------
    column_density : numpy.ndarray
        2D array of column densities
    """

    if abundance_grid is not None:
        density_grid = density_grid * abundance_grid
    
    column_density = np.sum(density_grid, axis=axis) * cell_size
    return column_density

def read_h5_data(filename, density_name='rho', abundance_name='abundance'):
    """
    Read density and abundance data from an HDF5 file.
    
    Parameters:
    -----------
    filename : str
        Path to the HDF5 file
    density_name : str
        Name of the density dataset
    abundance_name : str
        Name of the abundance dataset
    
    Returns:
    --------
    density_data : numpy.ndarray
        3D density array
    abundance_data : numpy.ndarray or None
        3D abundance array (None if not found)
    cell_size : float
        Calculated cell size
    attributes : dict
        Dictionary of attributes from the HDF5 file
    """
    with h5py.File(filename, 'r') as f:
        density_data = np.array(f[density_name])
        attributes = dict(f[density_name].attrs)
        cell_size = 3.0856e18 * np.max(np.array(f['x'])[:,0]) / len(np.array(np.array(f['x'])[:,0]))
        
        try:
            abundance_data = np.array(f[abundance_name])
        except KeyError:
            abundance_data = None
            print(f"Warning: Abundance dataset '{abundance_name}' not found - using density only")
            
    return density_data, abundance_data, cell_size, attributes

def save_results(output_file, column_density, axis_name, units='cm^-2', description=None):
    """
    Save column density data to an HDF5 file.
    
    Parameters:
    -----------
    output_file : str
        Path to the output HDF5 file
    column_density : numpy.ndarray
        2D array of column densities
    axis_name : str
        Name of the axis used for projection
    units : str
        Units of the column density
    description : str or None
        Description of the calculation
    """
    with h5py.File(output_file, 'w') as f:
        ds = f.create_dataset(f'column_density_{axis_name}', data=column_density)
        ds.attrs['units'] = units
        ds.attrs['projection_axis'] = axis_name
        if description:
            ds.attrs['description'] = description

def main():
    import argparse
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    

    # Set up command line arguments
    parser = argparse.ArgumentParser(description='Calculate column density from HDF5 file')
    parser.add_argument('input_file', help='Input HDF5 file containing density data')
    parser.add_argument('--rho', default='rho', help='Name of density dataset')
    parser.add_argument('--abundance', default='abundance', 
                       help='Name of abundance dataset (skip to use density only)')
    parser.add_argument('--output', help='Output HDF5 file prefix for results')
    parser.add_argument('--visualize', action='store_true', help='Show visualization')
    args = parser.parse_args()
    
    # Read data from HDF5 file
    density_data, abundance_data, cell_size, attrs = read_h5_data(
        args.input_file, 
        density_name=args.rho,
        abundance_name=args.abundance
    )
    print(f"Read density data with shape: {density_data.shape}")
    print(f"Calculated cell size: {cell_size:.2e} cm")
    if abundance_data is not None:
        print(f"Read abundance data with shape: {abundance_data.shape}")
    
    # Calculate column densities along all three axes
    results = {}
    for axis, axis_name in enumerate(['x', 'y', 'z']):
        col_dens = calculate_column_density(
            density_data, 
            abundance_grid=abundance_data,
            axis=axis, 
            cell_size=cell_size
        )
        results[axis_name] = col_dens
        print(f"Column density along {axis_name} axis shape: {col_dens.shape}")
        
        # Save results if output prefix specified
        if args.output:
            output_file = f"{args.output}_{axis_name}.h5"
            desc = f"Column density of {args.rho}"
            if abundance_data is not None:
                desc += f" * {args.abundance}"
            save_results(output_file, col_dens, axis_name, 
                        units=attrs.get('units', 'cm^-2'),
                        description=desc)
            print(f"Saved results to {output_file}")
    
    # Visualization
    if args.visualize:
        plt.figure(figsize=(20, 5))

        for i, (axis_name, col_dens) in enumerate(results.items()):
            plt.subplot(1, 3, i+1)
            plt.imshow(col_dens, cmap='gist_stern', origin='lower',norm=colors.LogNorm(vmin=1e15,vmax=1e20))
            plt.colorbar(label=f'Column Density ({attrs.get("units", r"cm$^{-2}$")})')
            title = f'Column Density along {axis_name.upper()}'
            if abundance_data is not None:
                title += f'\n({args.rho} × {args.abundance})'
            plt.title(title)
        
        plt.tight_layout()
        plt.savefig(f"{args.abundance}"+".png",bbox_inches='tight',dpi=300)
        plt.show()

if __name__ == "__main__":
    main()
