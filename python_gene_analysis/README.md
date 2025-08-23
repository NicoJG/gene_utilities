# GENE Analysis Utilities

A Python package for analyzing GENE (Gyrokinetic Electromagnetic Numerical Experiment) simulation data. Provides convenient functions for reading, visualizing, and diagnosing GENE output files.


## Main Features

- Read and parse GENE scan, geometry, and diagnostic files
- Visualize scan results and geometry overlays
- Run IDL-based diagnostics and extract results
- Interactive widgets for scan analysis


## Example Usage

Import functions directly from the package:

```python
from python_gene_analysis import (
    read_scan_log,
    read_ky_scan_data,
    plot_ky_scan,
    zprofile_iv,
    scan_analysis_widget
    ...
)

# Read scan data
scan_df = read_scan_log("/path/to/scan/")

# Plot ky scan
plot_ky_scan(["/path/to/scan1", "/path/to/scan2"])

# Run IDL diagnostic
zprofile_iv("/path/to/problem_dir", "/path/to/output_dir", run_id="run1")

# Use interactive widget (in Jupyter/IPython)
scan_analysis_widget(["/path/to/scan1", "/path/to/scan2"])
```
To use the IDL diagnostic automatically you should set `$GENEDIR` and the output of the IDL is always in a subdirectory of the scan directory called `idl_output`.

## Installation/Temporary Usage

You need to make python aware of this package. There are multiple ways to do this:d
1) Install this package into your Python environment:
```bash
pip install -e /path/to/python_gene_analysis
```
2) Add the parent directory to your shell environment variables either temporarily or for example in your `.bashrc`.
```bash
export PYTHONPATH="/path/to/parent:$PYTHONPATH"
```
3) Make it available in a python script (or notebook):
```python
import os
os.sys.path.append('/path/to/python_gene_analysis')
```

## Requirements
- numpy
- pandas
- matplotlib
- h5py (for cl_diag functionality)
- IDL installation (for automatic diagnostic generation)
