"""
GENE plotting - Collection of visualization functions for GENE simulation data.

Functions:
- plot_ky_scan()                   : ky scan plots with growth rates and frequencies
- plot_z_amplitude_on_geometry()   : z-amplitude overlaid on magnetic geometry
- plot_ballooning_potential()      : ballooning space potential vs poloidal angle
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from .gene_readers import read_ky_scan_data, read_z_amplitude_data, read_ballooning_data


def plot_ky_scan(scandirs, show=True, save_path=None, show_labels=False, figsize=(10,6)):
    """
    Create a ky scan plot showing growth rates and frequencies.
    
    Parameters:
    -----------
    scandirs : str, Path, or list of str/Path
        Single directory path or list of directory paths containing scan.log files
    show : bool, optional
        Whether to display the plot (default: True)
    save_path : str or Path, optional
        Path to save the plot (default: None)
    show_labels : bool, optional
        Whether to show (scan_nr, run) labels on data points (default: False)
    figsize : tuple, optional
        Figure size (default: (10, 6))
        
    Returns:
    --------
    tuple
        (fig, (ax1, ax2)) matplotlib figure and axes objects
    """
    # Read and prepare all data
    df_scan_log, params = read_ky_scan_data(scandirs)

    fig, ax1 = plt.subplots(figsize=figsize, layout='constrained')
    
    # Plot gamma on the left y-axis
    color1 = 'C0'
    ax1.set_xlabel(r'$k_y \rho_s$')
    ax1.set_ylabel('Growth Rate $\\gamma$', color=color1)
    line1 = ax1.plot(df_scan_log['kymin'], df_scan_log['gamma'], color=color1, marker='o', label='$\\gamma$')
    ax1.tick_params(axis='y', labelcolor=color1)
    ax1.grid(True, alpha=0.3)

    # Create second y-axis for omega
    ax2 = ax1.twinx()
    color2 = 'C1'
    ax2.set_ylabel('Frequency $\\omega$', color=color2)
    line2 = ax2.plot(df_scan_log['kymin'], df_scan_log['omega'], color=color2, marker='s', label='$\\omega$')
    ax2.tick_params(axis='y', labelcolor=color2)

    # Add tiny labels for each point showing (scan_nr, run) if enabled
    if show_labels and 'scan_nr' in df_scan_log.columns:
        for idx, row in df_scan_log.iterrows():
            label = f"({row['scan_nr']}, {row['run']})"
            # Add label to gamma points (on ax1)
            ax1.annotate(label, (row['kymin'], row['gamma']), ha='center', 
                        xytext=(5, 5), textcoords='offset points',
                        fontsize=6, alpha=0.7, color=color1)
            # Add label to omega points (on ax2)
            ax2.annotate(label, (row['kymin'], row['omega']), 
                        xytext=(5, -10), textcoords='offset points', ha='center',
                        fontsize=6, alpha=0.7, color=color2)

    # Combine legends
    lines = line1 + line2
    labels = [l.get_label() for l in lines]
    ax1.legend(lines, labels, loc='upper left')

    plt.title(f'$k_y \\rho_s$ scan for $a/L_n = {params["ions_omn"]}$ , $a/L_{{T_i}} = {params["ions_omt"]}$ , $a/L_{{T_e}} = {params["electrons_omt"]}$')
    if save_path:
        Path(save_path).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(save_path)
    if show:
        plt.show()
    return fig, (ax1, ax2)


def plot_z_amplitude_on_geometry(parameters_path, show=True, save_path=None, figsize=(13, 5), run_idl_automatically=False, geom_dir=None):
    """
    Plot z-amplitude overlaid on magnetic geometry.
    
    Parameters:
    -----------
    parameters_path : str or Path
        Path to the parameters file
    show : bool, optional
        Whether to display the plot (default: True)
    save_path : str or Path, optional
        Path to save the plot (default: None)
    figsize : tuple, optional
        Figure size (default: (13, 5))
    run_idl_automatically : bool, optional
        If True and zprofile file doesn't exist, automatically run IDL diagnostic
        
    Returns:
    --------
    tuple
        (fig, (ax1, ax2)) matplotlib figure and axes objects
    """
    # Read and prepare all data
    params, (kymin, gamma, omega), df_gist, df_all_vars, df_amplitude = read_z_amplitude_data(
        parameters_path, run_idl_automatically=run_idl_automatically, geom_dir=geom_dir
    )

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize, layout='constrained')

    # Plot 1: Amplitude and B(z)
    color1 = 'C1'
    ax1.set_xlabel('z')
    ax1.set_ylabel(r'$\sqrt{<|\phi|^2>}$', color=color1)
    line1 = ax1.plot(df_amplitude['z/qR'], df_amplitude['<|A|^2>^0.5'], color=color1, label=r'$\sqrt{<|\phi|^2>}$')
    ax1.tick_params(axis='y', labelcolor=color1)
    ax1.grid(True, alpha=0.3)

    ax1_twin = ax1.twinx()
    color2 = 'C0'
    ax1_twin.set_ylabel('$B$', color=color2)
    line2 = ax1_twin.plot(df_gist['z'], df_gist['B'], color=color2, label='B')
    ax1_twin.tick_params(axis='y', labelcolor=color2)

    ax1.set_title(r'Amplitude $\sqrt{<|\phi|^2>}$ and Magnetic Field $B$')

    # Combine legends for plot 1
    lines1 = line1 + line2
    labels1 = [l.get_label() for l in lines1]
    ax1.legend(lines1, labels1, loc='upper right')

    # Plot 2: Amplitude and K_x
    color3 = 'C1'
    ax2.set_xlabel('z')
    ax2.set_ylabel(r'$\sqrt{<|\phi|^2>}$', color=color3)
    line3 = ax2.plot(df_amplitude['z/qR'], df_amplitude['<|A|^2>^0.5'], color=color3, label=r'$\sqrt{<|\phi|^2>}$')
    ax2.tick_params(axis='y', labelcolor=color3)
    ax2.grid(True, alpha=0.3)

    ax2_twin = ax2.twinx()
    color4 = 'C2'
    ax2_twin.set_ylabel(r'$\kappa_\text{norm}$', color=color4)
    line4 = ax2_twin.plot(df_gist['z'], df_gist["k_norm"], color=color4, label=r'$\kappa_\text{norm}$')
    ax2_twin.tick_params(axis='y', labelcolor=color4)

    ax2.set_title(r'Amplitude $\sqrt{<|\phi|^2>}$ and Normal Curvature $\kappa_\text{norm}$')

    # Combine legends for plot 2
    lines2 = line3 + line4
    labels2 = [l.get_label() for l in lines2]
    ax2.legend(lines2, labels2, loc='upper right')

    fig.suptitle(f"$\\phi$ amplitude on magnetic geometry for\n$k_y \\rho_s = {kymin}$ , $\\gamma = {gamma}$ , $\\omega = {omega}$")

    if save_path:
        Path(save_path).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(save_path)
    if show:
        plt.show()
    return fig, (ax1, ax2)


def plot_ballooning_potential(parameters_path, show=True, save_path=None, figsize=(10, 6), run_idl_automatically=False):
    """
    Plot the ballooning space potential φ vs poloidal angle.
    
    Parameters:
    -----------
    parameters_path : str or Path
        Path to the parameters file
    show : bool, optional
        Whether to display the plot (default: True)
    save_path : str or Path, optional
        Path to save the plot (default: None)
    figsize : tuple, optional
        Figure size (default: (10, 6))
    run_idl_automatically : bool, optional
        If True and file doesn't exist, automatically run IDL diagnostic to generate it
        
    Returns:
    --------
    tuple
        (fig, ax) matplotlib figure and axis objects
    """
    # Read and prepare all data
    df_ball, params, (kymin, gamma, omega) = read_ballooning_data(
        parameters_path, run_idl_automatically=run_idl_automatically
    )
    
    # Create the plot
    fig, axs = plt.subplots(2, 1, figsize=figsize, layout='constrained')

    plt.sca(axs[0])

    plt.plot(df_ball['theta'], df_ball['Re_phi'], color='tab:red', 
                linewidth=1.5, label=r'$\Re(\phi)$')
    plt.plot(df_ball['theta'], df_ball['Im_phi'], color='tab:blue', 
                linewidth=1.5, label=r'$\Im(\phi)$')
    plt.plot(df_ball['theta'], df_ball['abs_phi'], 'k-', 
            linewidth=2, label=r'$|\phi|$')

    plt.xlabel(r'$\theta / \pi$')
    plt.ylabel(r'$\phi$')
    plt.title(f'Ballooning space potential $\\phi$ for\n$k_y = {kymin}$ , $\\gamma = {gamma}$ , $\\omega = {omega}$')
    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.sca(axs[1])

    plt.plot(df_ball['theta'], np.abs(df_ball['Re_phi']), color='tab:red', 
                linewidth=1.5, label=r'$\Re(\phi)$')
    plt.plot(df_ball['theta'], np.abs(df_ball['Im_phi']), color='tab:blue', 
                linewidth=1.5, label=r'$\Im(\phi)$')
    plt.plot(df_ball['theta'], np.abs(df_ball['abs_phi']), 'k-', 
            linewidth=2, label=r'$|\phi|$')

    plt.xlabel(r'$\theta / \pi$')
    plt.ylabel(r'$\phi$')
    plt.yscale("log")


    if save_path:
        Path(save_path).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(save_path)
    if show:
        plt.show()
        
    return fig, axs
