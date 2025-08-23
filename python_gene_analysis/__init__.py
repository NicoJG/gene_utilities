"""
GENE Analysis Utilities

A collection of utility functions for analyzing GENE simulation data.

Modules:
- gene_readers: Functions for reading GENE data files
- gene_plotting: Functions for plotting GENE data
- cl_diag: IDL diagnostic interface utilities
- gene_analysis_utils: Interactive widgets for scan analysis
"""

# Import all main functions for convenience
from .gene_readers import (
    read_scan_log,
    read_ky_scan_data,
    read_z_amplitude_data,
    read_ballooning_data,
    read_gist_file,
    read_zprofile,
    read_parameter_file
)

from .gene_plotting import (
    plot_ky_scan,
    plot_z_amplitude_on_geometry,
    plot_ballooning_potential
)

# Interactive widgets
from .interactive_widgets import scan_analysis_widget

from .cl_diag import (
    zprofile_iv,
    ball_iv,
    ball_ev,
    timetrace,
    run_diagnostic_iv,
    write_cl_diag
)

__all__ = [
    # gene_readers functions
    'read_scan_log',
    'read_ky_scan_data', 
    'read_z_amplitude_data',
    'read_ballooning_data',
    'read_gist_file',
    'read_zprofile',
    'read_parameter_file',
    # gene_plotting functions
    'plot_ky_scan',
    'plot_z_amplitude_on_geometry',
    'plot_ballooning_potential',
    # cl_diag functions
    'zprofile_iv',
    'ball_iv',
    'ball_ev',
    'timetrace',
    'run_diagnostic_iv',
    'write_cl_diag'
    # interactive widgets
    'scan_analysis_widget'
]