"""
GENE file readers - Collection of utility functions for reading GENE data files.

Functions:
- read_scan_log()          : Read scan.log files and extract eigenvalues
- read_ky_scan_data()      : Complete ky scan data with parameters  
- read_z_amplitude_data()  : Data for z-amplitude geometry plots
- read_ballooning_data()   : Data for ballooning potential plots
- read_gist_file()         : Read GIST geometry files
- read_zprofile()          : Read zprofile diagnostic files
- read_parameter_file()    : Parse GENE parameter files
"""

import numpy as np
import pandas as pd
from pathlib import Path
import re
from . import cl_diag


def read_scan_log(scandirs, extract_scan_numbers=True):
    """
    Read and process GENE scan.log file(s) to extract eigenvalues.
    
    This function reads scan.log files from one or multiple directories,
    parses the eigenvalue data, and splits it into gamma (growth rate) 
    and omega (frequency) components. It dynamically adapts to the 
    variables being scanned over by parsing the header line.
    
    Parameters:
    -----------
    scandirs : str, Path, or list of str/Path
        Single directory path or list of directory paths containing scan.log files
    extract_scan_numbers : bool, optional
        If True, attempts to extract scan numbers from directory paths (default: True)
        
    Returns:
    --------
    pd.DataFrame
        DataFrame with columns: run, [scan_variables...], gamma, omega, path, [scan_nr]
        The scan variables are automatically detected from the header line.
        If multiple directories are provided, data is concatenated with 
        reset index, excluding the first row of subsequent files.
        If extract_scan_numbers=True, includes a 'scan_nr' column.
    """
    def _extract_scan_nr(path_str):
        """Extract scan number from path string."""
        match = re.search(r'scanfiles(\d+)', str(path_str))
        if match:
            return int(match.group(1))
        else:
            # Try other common patterns
            match = re.search(r'scan[_-]?(\d+)', str(path_str), re.IGNORECASE)
            if match:
                return int(match.group(1))
            else:
                return str(path_str)  # Use full path if pattern doesn't match
    
    # Convert to Path object if string
    if isinstance(scandirs, (str, Path)):
        scandirs = Path(scandirs)
        
    # Handle single directory case
    if isinstance(scandirs, Path):
        # Read the header line to extract variable names
        with open(scandirs / "scan.log", 'r') as f:
            header_line = f.readline().strip()
        
        # Extract variable names using regex
        # Pattern matches: variable_name followed by whitespace and numbers, before /Eigenvalue
        pattern = r'\|\s*(\w+)\s+[\d\.\-e\+]+(?=.*?/Eigenvalue|\s*$)'
        variable_names = re.findall(pattern, header_line)
        
        # Build column names
        column_names = ["run"] + variable_names + ["eigenvalues"]
        
        # Read and process data
        df = pd.read_csv(scandirs / "scan.log", sep="|", skiprows=1, header=None, names=column_names)
        
        # Split eigenvalues into gamma and omega using pandas vectorized operations
        eigenvalue_df = df['eigenvalues'].str.strip().str.split(expand=True).astype(float)
        df = df.assign(
            gamma=eigenvalue_df[0],
            omega=eigenvalue_df[1], 
            path=str(scandirs.name)
        ).drop('eigenvalues', axis=1)

        # Add scan number if requested
        if extract_scan_numbers:
            df['scan_nr'] = _extract_scan_nr(scandirs.name)

        # Sort by scan variables (excluding run, gamma, omega, path, and optionally scan_nr)
        sort_cols = df.columns[1:-3].tolist()  # Exclude run, gamma, omega, path
        if extract_scan_numbers:
            sort_cols = sort_cols[:-1]  # Also exclude scan_nr
        df.sort_values(by=sort_cols, inplace=True)

        return df
    
    # Handle multiple directories case using recursion
    elif isinstance(scandirs, list):
        if not scandirs:
            raise ValueError("scandirs list cannot be empty")
            
        # Process all directories
        dfs = []
        for i, scandir in enumerate(scandirs):
            dfs.append(read_scan_log(scandir, extract_scan_numbers=extract_scan_numbers))
            
        # Concatenate all dataframes
        result_df = pd.concat(dfs, axis=0, ignore_index=True)
        
        # Sort by scan variables if we have them
        if len(result_df) > 0:
            sort_cols = result_df.columns[1:-3].tolist()  # Exclude run, gamma, omega, path
            if extract_scan_numbers and 'scan_nr' in result_df.columns:
                sort_cols = sort_cols[:-1]  # Also exclude scan_nr
            if sort_cols:  # Only sort if we have scan variables
                result_df.sort_values(by=sort_cols, inplace=True)
        
        return result_df
    
    else:
        raise TypeError(f"scandirs must be str, Path, or list of str/Path, got {type(scandirs)}")


def read_ky_scan_data(scandirs):
    """
    Read and prepare complete ky scan data including scan log and parameters.
    
    This is a higher-level function that combines scan log reading with parameter
    extraction, providing all data needed for ky scan plotting.
    
    Parameters:
    -----------
    scandirs : str, Path, or list of str/Path
        Single directory path or list of directory paths containing scan.log files
        
    Returns:
    --------
    tuple[pd.DataFrame, dict]
        - DataFrame with scan data (sorted by kymin, includes scan_nr)
        - Dictionary with parameters from the first scan directory
    """
    # Read scan log data with scan number extraction
    df_scan_log = read_scan_log(scandirs, extract_scan_numbers=True)
    
    # Sort by kymin for plotting
    if 'kymin' in df_scan_log.columns:
        df_scan_log.sort_values(by="kymin", inplace=True)
    
    # Get parameters from the first scan directory
    if isinstance(scandirs, list):
        first_dir = Path(scandirs[0])
    else:
        first_dir = Path(scandirs)
    
    params = read_parameter_file(first_dir / "parameters")
    
    return df_scan_log, params


def read_z_amplitude_data(parameters_path, run_idl_automatically=False, geom_dir=None):
    """
    Read and prepare data for z-amplitude plotting on geometry.
    
    This function reads parameters, omega data, GIST geometry, and zprofile data
    needed for creating z-amplitude plots overlaid on magnetic geometry.
    
    Parameters:
    -----------
    parameters_path : str or Path
        Path to the parameters file
    run_idl_automatically : bool, optional
        If True and zprofile file doesn't exist, automatically run IDL diagnostic
        
    Returns:
    --------
    tuple[dict, tuple, pd.DataFrame, pd.DataFrame, pd.DataFrame]
        - params: Parameter dictionary
        - (kymin, gamma, omega): Eigenvalue data
        - df_gist: GIST geometry DataFrame  
        - df_all_vars: Complete zprofile DataFrame
        - df_amplitude: Processed amplitude DataFrame for plotting
    """
    parameters_path = Path(parameters_path)
    
    # Read run parameters
    params = read_parameter_file(parameters_path)
    
    # Read omega data
    runnumber = parameters_path.name.replace('parameters_', '')
    kymin, gamma, omega = np.genfromtxt(parameters_path.parent / f"omega_{runnumber}")
    
    # Read GIST geometry data
    if geom_dir is None:
        gist_path = Path(params["geomdir"]) / params["geomfile"]
    else:
        gist_path = Path(geom_dir) / params["geomfile"]
    df_gist, dict_gist = read_gist_file(gist_path)
    
    # Read zprofile data with case variations
    zprofile_lower = parameters_path.parent / f"idl_output/zprofileions_{runnumber}.dat"
    zprofile_upper = parameters_path.parent / f"idl_output/ZPROFILEions_{runnumber}.dat"
    
    if zprofile_lower.exists():
        zprofile_file = zprofile_lower
    elif zprofile_upper.exists():
        zprofile_file = zprofile_upper
    else:
        # Try to automatically generate if requested
        if run_idl_automatically:
            zprofile_file = zprofile_lower  # Will be generated by read_zprofile
        else:
            raise FileNotFoundError(f"No zprofile ions file found for run {runnumber}. Checked:\n"
                                  f"  {zprofile_lower}\n  {zprofile_upper}\n"
                                  f"Set run_idl_automatically=True to generate automatically.")
    
    # Read zprofile data (single call)
    df_all_vars = read_zprofile(zprofile_file, run_idl_automatically=run_idl_automatically, parameters_path=parameters_path)
    
    # Extract phi amplitude for plotting
    df_amplitude = pd.DataFrame({
        'z/qR': df_all_vars['z/qR'],
        '<|A|^2>^0.5': df_all_vars['<|phi|^2>^1/2']
    })
    
    return params, (kymin, gamma, omega), df_gist, df_all_vars, df_amplitude


def read_ballooning_data(parameters_path, run_idl_automatically=False):
    """
    Read and prepare data for ballooning potential plotting.
    
    This function reads the ballooning data file and parameters needed
    for creating ballooning potential plots.
    
    Parameters:
    -----------
    parameters_path : str or Path
        Path to the parameters file
    run_idl_automatically : bool, optional
        If True and ballooning file doesn't exist, automatically run IDL diagnostic
        
    Returns:
    --------
    tuple[pd.DataFrame, dict, tuple]
        - df_ball: Ballooning data DataFrame with calculated abs_phi
        - params: Parameter dictionary  
        - (kymin, gamma, omega): Eigenvalue data
    """
    parameters_path = Path(parameters_path)
    
    # Extract run number from parameters file
    runnumber = parameters_path.name.replace('parameters_', '')
    
    # Read ballooning data with case variations
    ball_upper = parameters_path.parent / f"idl_output/BALLions_{runnumber}.dat"
    ball_lower = parameters_path.parent / f"idl_output/ballions_{runnumber}.dat"
    
    if ball_upper.exists():
        ball_file = ball_upper
    elif ball_lower.exists():
        ball_file = ball_lower
    else:
        # Try to automatically generate if requested
        if run_idl_automatically:
            # Create output directory if it doesn't exist
            output_dir = parameters_path.parent / "idl_output"
            output_dir.mkdir(exist_ok=True)
            
            cl_diag.ball_iv(prob_dir=str(parameters_path.parent), 
                           out_dir=str(output_dir),
                           run_id=runnumber)
            
            # After generation, assume lowercase (like zprofile)
            ball_file = ball_lower
        else:
            raise FileNotFoundError(f"No ballooning data file found for run {runnumber}. Checked:\n"
                                  f"  {ball_upper}\n  {ball_lower}\n"
                                  f"Set run_idl_automatically=True to generate automatically.")
    
    # Read the data with automatic numeric conversion (single call)
    df_ball = pd.read_csv(ball_file,
                          sep=r'\s+', 
                          skiprows=10,
                          header=None,
                          names=['theta', 'Re_phi', 'Im_phi', 'theta_pol'],
                          dtype=float,
                          encoding='iso-8859-1')
    
    # Calculate absolute value using complex representation
    df_ball['abs_phi'] = np.abs(df_ball['Re_phi'] + 1j * df_ball['Im_phi'])
    
    # Read parameters and omega data
    params = read_parameter_file(parameters_path)
    kymin, gamma, omega = np.genfromtxt(parameters_path.parent / f"omega_{runnumber}")
    
    return df_ball, params, (kymin, gamma, omega)


def read_gist_file(filepath):
    """
    Read a GIST geometry file and extract both data and parameters.
    
    GIST files contain a Fortran namelist section with parameters followed by
    numerical data in columns. This function parses both sections and validates
    that the file was created through GeneTools.jl.
    
    Parameters:
    -----------
    filepath : str or Path
        Path to the GIST file
        
    Returns:
    --------
    tuple[pd.DataFrame, dict]
        - DataFrame with the numerical data (9 columns expected)
        - Dictionary with the parameters from the namelist section
    
    Raises:
    -------
    ValueError
        If the file doesn't have exactly 9 columns (not created through GeneTools.jl)
    """
    def _parse_value(value_str):
        """Convert string value to appropriate numeric type."""
        try:
            return float(value_str) if ('.' in value_str or 'e' in value_str.lower()) else int(value_str)
        except ValueError:
            return value_str
    
    def _extract_parameters(line):
        """Extract key-value pairs from a parameter line."""
        if '=' not in line:
            return {}
        
        key, value = line.split('=', 1)
        key, value = key.strip(), value.strip()
        
        # Extract numerical values
        numbers = re.findall(r'-?\d+\.?\d*(?:[eE][+-]?\d+)?', value)
        
        # Handle comma-separated keys (e.g., "s0, alpha0")
        if ',' in key:
            keys = [k.strip() for k in key.split(',')]
            if len(keys) == len(numbers):
                return {k: _parse_value(num) for k, num in zip(keys, numbers)}
            else:
                return {key: value}  # Fallback if mismatch
        
        # Handle single key
        if len(numbers) == 1:
            return {key: _parse_value(numbers[0])}
        elif len(numbers) > 1:
            return {key: [_parse_value(x) for x in numbers]}
        else:
            return {key: value}
    
    with open(Path(filepath), 'r') as f:
        lines = f.readlines()
    
    # Find namelist end and parse parameters
    namelist_end = next(i for i, line in enumerate(lines) if line.strip() == '/')
    parameters = {}
    coordinate_system = None
    
    for line in lines[1:namelist_end]:
        line = line.strip()
        
        # Detect coordinate system
        if line.startswith('!') and ('PEST' in line.upper() or 'BOOZER' in line.upper()):
            coordinate_system = 'PEST' if 'PEST' in line.upper() else 'BOOZER'
        
        # Parse comment parameters (!key = value)
        elif line.startswith('!') and '=' in line:
            parameters.update(_extract_parameters(line[1:].strip()))
        
        # Parse regular parameters (key = value)
        elif '=' in line and not line.startswith('!'):
            clean_line = line.split('!')[0].strip()  # Remove inline comments
            if clean_line:
                parameters.update(_extract_parameters(clean_line))
    
    if coordinate_system:
        parameters['coordinate_system'] = coordinate_system
    
    # Parse numerical data
    data_rows = [[float(x) for x in line.split()] 
                 for line in lines[namelist_end + 1:] if line.strip()]
    
    # Validate GeneTools.jl format (9 columns)
    if data_rows and len(data_rows[0]) != 9:
        raise ValueError(f"GIST file has {len(data_rows[0])} columns but GeneTools.jl creates "
                        f"files with exactly 9 columns. This file was likely not created through GeneTools.jl")
    
    # Create DataFrame with proper column names
    columns = ['g_xx', 'g_xy', 'g_yy', 'B', 'jacobian', 'K_2', 'K_1', 'dBdz', 'last_col']
    df = pd.DataFrame(data_rows, columns=columns)
    
    # Apply corrections and calculate derived quantities
    df["K_1"] = -df["K_1"]  # GIST stores negative of K1
    df["K_2"] = df["K_2"] - parameters["my_dpdx"]/(2*df["B"])

    # add z coordinate as first column
    df.insert(0, 'z', np.linspace(-np.pi*parameters["n_pol"],np.pi*parameters["n_pol"], parameters["gridpoints"], endpoint=False))

    # Calculate normal and geodesic curvature
    df["k_geo"] = df["K_1"]/df["g_xx"]**0.5
    df["k_norm"] = (df["g_xx"]*df["K_2"] - df["g_xy"]*df["k_geo"])/(df["B"]*df["g_xx"]**0.5)

    return df, parameters


def read_zprofile(filepath, run_idl_automatically=False, parameters_path=None):
    """
    Read zprofile file and merge all variables into one DataFrame.
    
    Parameters:
    -----------
    filepath : str or Path
        Path to the zprofile file to read
    run_idl_automatically : bool, optional
        If True and file doesn't exist, automatically run IDL diagnostic to generate it
    parameters_path : str or Path, optional
        Path to parameters file, needed for automatic IDL execution
    
    Returns:
    --------
    pd.DataFrame: Single DataFrame with z/qR column and columns for each variable
    """
    
    # Check if file exists, if not and auto-run is enabled, generate it
    filepath = Path(filepath)
    if not filepath.exists() and run_idl_automatically:
        if parameters_path is None:
            raise ValueError("parameters_path must be provided when run_idl_automatically=True")
            
        parameters_path = Path(parameters_path)
        run_id = parameters_path.name.replace("parameters_", "")
        
        # Create output directory if it doesn't exist
        output_dir = parameters_path.parent / "idl_output"
        output_dir.mkdir(exist_ok=True)
        
        cl_diag.zprofile_iv(prob_dir=str(parameters_path.parent), 
                           out_dir=str(output_dir),
                           run_id=run_id)
    
    with open(filepath, 'r', encoding='iso-8859-1') as f:
        lines = f.readlines()
    
    blocks = {}
    current_var = None
    current_data = []
    
    for line in lines:
        if line.startswith('# variable A ='):
            # Save previous variable
            if current_var and current_data:
                blocks[current_var] = current_data
            # Start new variable
            current_var = line.split('=')[1].strip()
            current_data = []
        elif not line.startswith('#') and line.strip():
            try:
                values = line.split()
                if len(values) == 4:
                    current_data.append([float(val) for val in values])
            except ValueError:
                pass
    
    # Don't forget the last variable
    if current_var and current_data:
        blocks[current_var] = current_data
    
    if not blocks:
        return pd.DataFrame()
    
    # Use first variable for z-coordinates
    first_data = list(blocks.values())[0]
    df = pd.DataFrame({'z/qR': [row[0] for row in first_data]})
    
    # Add each variable with mathematical notation
    for var_name, data in blocks.items():
        df[f'<|{var_name}|>'] = [row[1] for row in data]
        df[f'<|{var_name}|^2>^1/2'] = [row[2] for row in data]
        df[f'σ_{var_name}'] = [row[3] for row in data]
    
    return df


def read_parameter_file(filepath):
    """
    Read a GENE parameter file into a Python dictionary.
    Ignores namelist tags (&tag) and parses variable = value ! comment lines.
    For &species namelists, prefixes variables with species name (e.g., ions_omn, electrons_omt).
    
    Parameters:
    -----------
    filepath : str or Path
        Path to the parameter file
        
    Returns:
    --------
    dict : Dictionary with variable names as keys and values as parsed values
    """
    params = {}
    current_species = None
    in_species_block = False
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            
            # Skip empty lines and comment-only lines
            if not line or line.startswith('!'):
                continue
                
            # Check for namelist start
            if line.startswith('&'):
                if line.startswith('&species'):
                    in_species_block = True
                    current_species = None  # Reset species name
                else:
                    in_species_block = False
                    current_species = None
                continue
                
            # Check for namelist end
            if line.startswith('/'):
                in_species_block = False
                current_species = None
                continue
                
            # Look for variable = value pattern
            if '=' in line:
                # Split on first '=' to handle cases where value contains '='
                var_part, value_part = line.split('=', 1)
                var_name = var_part.strip()
                
                # Remove comments (everything after '!')
                if '!' in value_part:
                    value_part = value_part.split('!')[0]
                
                value_str = value_part.strip()
                
                # Parse the value - try different types
                try:
                    # Try integer first
                    if '.' not in value_str and 'E' not in value_str.upper():
                        value = int(value_str)
                    else:
                        # Try float
                        value = float(value_str)
                except ValueError:
                    # Handle boolean values
                    if value_str.upper() in ['T', 'TRUE', '.TRUE.']:
                        value = True
                    elif value_str.upper() in ['F', 'FALSE', '.FALSE.']:
                        value = False
                    else:
                        # Keep as string, remove quotes if present
                        value = value_str.strip('\'"')
                
                # Handle species block variables
                if in_species_block:
                    if var_name == 'name':
                        # Store the species name for prefixing other variables
                        current_species = value
                        params[var_name] = value  # Also store the name itself
                    else:
                        # Prefix variable with species name if we have one
                        if current_species:
                            prefixed_name = f"{current_species}_{var_name}"
                            params[prefixed_name] = value
                        else:
                            # Fallback if name wasn't found yet
                            params[var_name] = value
                else:
                    # Regular variable outside species block
                    params[var_name] = value
    
    return params
