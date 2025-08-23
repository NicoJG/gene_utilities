import ipywidgets as widgets
from IPython.display import display
import os
import glob
import re
from pathlib import Path
from .gene_plotting import plot_z_amplitude_on_geometry, plot_ballooning_potential

def scan_analysis_widget(scandir_paths, run_idl_automatically=False, geom_dir=None):
    """
    Create an interactive widget to select runs and plot diagnostics for scanfiles.

    Features:
    ---------
    - Dropdown to select scan folder
    - Dropdown to select run
    - Tick boxes to select which diagnostics to plot (Ballooning, Phi Z Amplitude, etc.)
    - Easy to add more diagnostics by editing the diagnostic_options list

    Parameters:
    -----------
    scandir_paths : Path or list of Path
        Path to a scanfiles directory, or list of paths to multiple scanfiles directories
    """

    # Handle both single path and list of paths
    if not isinstance(scandir_paths, list):
        scandir_paths = [scandir_paths]

    # Build dictionary of all available data
    all_data = {}  # {folder_name: {run_num: file_path}}

    for scandir_path in scandir_paths:
        scandir_path = Path(scandir_path)
        folder_name = scandir_path.name

        # Find parameter files with pattern parameters_XXXX (digits only)
        parameter_files = sorted(scandir_path.glob("parameters_*"))
        pattern = re.compile(r'^parameters_(\d+)$')

        run_files = {}
        for pfile in parameter_files:
            filename = pfile.name
            match = pattern.match(filename)
            if match:
                run_num = match.group(1)
                run_files[run_num] = pfile

        if run_files:
            all_data[folder_name] = run_files
    
    if not all_data:
        print(f"No parameter files found in any of the provided directories")
        return
    
    # Create widgets
    folder_names = list(all_data.keys())
    
    # Always show folder dropdown, even if only one folder
    folder_dropdown = widgets.Dropdown(
        options=folder_names,
        value=folder_names[0],
        description='Folder:'
    )
    
    # Initialize run dropdown with first folder's runs
    initial_runs = list(all_data[folder_names[0]].keys())
    run_dropdown = widgets.Dropdown(
        options=initial_runs,
        value=initial_runs[0] if initial_runs else None,
        description='Run:'
    )
    
    # Diagnostic selection checkboxes
    diagnostic_options = [
        ('Ballooning', 'ballooning'),
        ('Phi Z Amplitude', 'phi_z_amplitude'),
        # Add more diagnostics here as needed
    ]
    diagnostic_checkboxes = {key: widgets.Checkbox(value=(key=='phi_z_amplitude'), description=label)
                             for label, key in diagnostic_options}

    diagnostics_box = widgets.VBox(list(diagnostic_checkboxes.values()))

    output = widgets.Output()
    
    def update_run_options(change):
        """Update run dropdown when folder selection changes"""
        selected_folder = change['new']
        available_runs = list(all_data[selected_folder].keys())
        run_dropdown.options = available_runs
        if available_runs:
            run_dropdown.value = available_runs[0]
    
    def plot_selected_run(change=None):
        """Callback for run dropdown or diagnostic selection changes"""
        with output:
            output.clear_output(wait=True)
            selected_folder = folder_dropdown.value
            selected_run = run_dropdown.value
            print(f"Plotting {selected_folder}/run {selected_run}...")
            try:
                if diagnostic_checkboxes['phi_z_amplitude'].value:
                    plot_z_amplitude_on_geometry(all_data[selected_folder][selected_run], show=True, run_idl_automatically=run_idl_automatically, geom_dir=geom_dir)
                if diagnostic_checkboxes['ballooning'].value:
                    plot_ballooning_potential(all_data[selected_folder][selected_run], show=True, run_idl_automatically=run_idl_automatically)
                # Add more diagnostics here as needed
            except Exception as e:
                print(f"Error plotting {selected_folder}/run {selected_run}: {e}")
    
    # Setup observers
    folder_dropdown.observe(update_run_options, names='value')
    run_dropdown.observe(plot_selected_run, names='value')
    for cb in diagnostic_checkboxes.values():
        cb.observe(plot_selected_run, names='value')
    
    # Display widgets (always show folder dropdown)
    display(folder_dropdown, run_dropdown, diagnostics_box, output)
    
    # Plot initial selection
    with output:
        plot_selected_run()