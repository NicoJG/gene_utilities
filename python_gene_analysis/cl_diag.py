"""
GENE IDL diagnostics - Command-line interface to GENE IDL diagnostic tools.

Functions:
- write_cl_diag()        : Write command line diagnostic parameters file
- run_diagnostic_iv()    : Run any IDL diagnostic with specified parameters
- ball_iv()              : Run ballooning diagnostic via IDL
- zprofile_iv()          : Run zprofile diagnostic via IDL  
- ball_ev()              : Run eigenspectrum ballooning diagnostic
- timetrace()            : Run timetrace diagnostic via IDL

Originally written by Michael Gerard (https://gitlab.com/mjgerard2/my_gene)
"""

import os
import shlex
import subprocess
import h5py as hf
import numpy as np


def get_diag_dir():
    """Get the diagnostics directory from $GENEDIR, or raise an error if not set."""
    genedir = os.environ.get('GENEDIR')
    if genedir is not None:
        return os.path.join(genedir, 'diagnostics')
    raise EnvironmentError(
        "$GENEDIR is not set. Please set the GENEDIR environment variable, e.g. 'export GENEDIR=/path/to/gene' in your shell. "
        "Alternatively, set it in Python before importing this module: "
        "import os; os.environ['GENEDIR'] = '/path/to/gene'"
    )


def write_cl_diag(cl_dict=None):
    """ Write the command line diagnostic parameters file.

    Parameters
    ----------
    cl_dict: dict
        Dictionary of command line parameters.
    """
    diag_dir = get_diag_dir()
    cl_params_file = os.path.join(diag_dir, 'cl_params')
    with open(cl_params_file, 'w', encoding='utf-8') as file:
        file.write('{}/ ; '.format(cl_dict['data path']) +
                   'data path, slash at end mandatory\n')
        file.write('{}/ ; '.format(cl_dict['output path']) +
                   'output path, slash at end mandatory\n')
        file.write('{} ; start time\n'.format(cl_dict['start time']))
        file.write('{} ; end time\n'.format(cl_dict['end time']))
        file.write('{} ; sparse factor\n'.format(cl_dict['sparse factor']))
        file.write('{} ; '.format(cl_dict['run ID']) +
                   'run number(s) or label(s)\n')
        file.write('{} ; '.format(cl_dict['species']) +
                   'species number (same order as in parameters file)\n')
        file.write('{} ; '.format(cl_dict['diag']) +
                   'diagnostic name (can specify only one, only mom)\n')
        file.write('---------------------------------------------------\n')
        file.write('Do not change the number of lines above this point!\n')
        file.write('Comments may be modified/extended/removed\n')
        file.write('---------------------------------------------------\n')


def run_diagnostic_iv(prob_dir, out_dir, run_id, diag_type, times=[0, 1e4], species=0, sparse_factor=1, verbose=True):
    """
    Run any diagnostic through the idl command-line tools.
    It is typically helpful to set 'gui.out.only_last = 1' around line 40 in /gene/diagnostics/cl_diag.pro.
    Also, for ASCII outputs, change 'postscript' to 'data: ASCII' when defining the ps_format variable in the same file.

    Parameters
    ----------
    prob_dir: str
        Global directory to the output problem directory.
    out_dir: str
        Global directory to where the data is output.
    run_id: int
        Run ID in the problem directory.
    diag_type: str
        Type of diagnostic to run (e.g., 'zprofile', 'ball', 'spectra', etc.)
    times: tuple, optional
        The start and end times for the diag parameters. This can be ignored
        if gui.out.only_last is set to 1. Default is [0, 1e4].
    species: int, optional
        Species index for the diagnostic. Default is 0.
    sparse_factor: int, optional
        Sparse factor for the diagnostic. Default is 1.
    verbose: bool, optional
        Whether to print IDL output for debugging. Default is True.

    Returns
    -------
    tuple
        (stdout, stderr) from the IDL process for debugging
    """
    # write command line parameters file #
    cl_dict = {'data path': prob_dir,
               'output path': out_dir,
               'start time': times[0],
               'end time': times[-1],
               'sparse factor': sparse_factor,
               'run ID': run_id,
               'species': species,
               'diag': diag_type}
    write_cl_diag(cl_dict)

    # execute IDL diagnostic in child process #
    diag_dir = get_diag_dir()
    wdir = os.getcwd()
    os.chdir(diag_dir)
    proc_idl = subprocess.Popen(shlex.split('idl -quiet cl_diag.pro'),
                                stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                text=True)
    proc_out, proc_err = proc_idl.communicate(input='exit')
    os.chdir(wdir)
    if proc_err:
        raise RuntimeError(f"IDL call failed for diagnostic '{diag_type}'.\nSTDERR:\n{proc_err}")
    if verbose:
        print(f"IDL {diag_type} diagnostic stdout:")
        print(proc_out)
    return proc_out, proc_err


def ball_iv(prob_dir, out_dir, run_id, times=[0, 1e4]):
    return run_diagnostic_iv(prob_dir, out_dir, run_id, 'ball', times)


def zprofile_iv(prob_dir, out_dir, run_id, times=[0, 1e4]):
    return run_diagnostic_iv(prob_dir, out_dir, run_id, 'zprofile', times)


def ball_ev(prob_dir, out_dir, run_id, n_ev):
    """
    Run the ballooning diagnostic through the idl command-line tools for an eigenspectrum computation.
    If 'gui.out.only_last = 1' is defined anywhere in /gene/diagnostics/cl_diag.pro, then you should comment that line out.
    Also, ASCII outputs are required for reading the individual files that are then converted to HDF5.
    Therefore, change 'postscript' to 'data: ASCII' when defining the ps_format variable in the same file.

    Parameters
    ----------
    prob_dir: str
        Global directory to the output problem directory.
    out_dir: str
        Global directory to where the data is output.
    run_id: int
        Run ID in the problem directory.
    n_ev: int
        The number of eigenvectors.
    """
    diag_dir = get_diag_dir()
    # define global path to ASCII formatted output file #
    out_file = os.path.join(out_dir, f'ballions_{run_id}.dat')
    # write command line parameters file #
    cl_dict = {'data path': prob_dir,
               'output path': out_dir,
               'sparse factor': 1,
               'run ID': run_id,
               'species': 0,
               'diag': 'ball'}

    # change to diagnostic directory #
    wdir = os.getcwd()
    os.chdir(diag_dir)
    # loop over eigenvectors #
    for i in range(1, n_ev+1):
        # write command line diagnostic parameters #
        cl_dict['start time'] = i
        cl_dict['end time'] = i
        write_cl_diag(cl_dict)
        # execute cl_diag.pro script #
        proc_idl = subprocess.Popen(shlex.split('idl cl_diag.pro'),
                                    stdin=subprocess.PIPE,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    text=True)
        proc_out, proc_err = proc_idl.communicate(input='exit')
        if proc_err:
            raise RuntimeError(f"IDL call failed in ball_ev (eigenvector {i}).\nSTDERR:\n{proc_err}")
        # read out ASCII formatted output into numpy array #
        with open(out_file, 'r') as f:
            lines = f.readlines()
            for l, line in enumerate(lines):
                line_strip = line.strip().split()
                if line_strip[0] == 'theta_M':
                    l += 1
                    break
            if i == 1:
                ball_data = np.empty((n_ev, len(lines[l::]), 4))
            for l, line in enumerate(lines[l::]):
                ball_data[i-1, l] = [float(x) for x in
                                          line.strip().split()]
        #os.remove(out_file)
    # save eigenvector data in HDF5 format #
    #save_file = os.path.join(out_dir, f'ballions_{run_id}.h5')
    #with hf.File(save_file, 'w') as hf_:
    #    hf_.create_dataset('ball data', data=ball_data)


def timetrace(dirs, run_id, times, var_idx, kx_idx, ky_idx):
    """ Run the timetrace diagnostic through the idl command-line tools. For
    ASCII outputs, change \'postscript\' to \'data: ASCII\' when defining the
    ps_format variable in the same file.

    Parameters
    ----------
    dirs: dict
        Contains the following global directories:
            problem: directory where the GENE data is stored
            output: directory where the post-processing data will be output
            diagnostic: GENE diagnostic directory (optional, ignored, always uses GENEDIR)
    run_id: int or list
        Run ID in the problem directory.
    times: tuple
        The start and end times for the diag parameters.
    var_idx: list
        Indices of the variables to be analyzed.
    kx_idx: list
        Indices of the radial wavenumbers to be output.
    ky_idx: list
        Indices of the binormal wavenumbers to be output.
    """
    # number of fluctuation arrays #
    nvar = len(var_idx)
    nkx = len(kx_idx)
    nky = len(ky_idx)

    # check if run_id is an integer or a list #
    if isinstance(run_id, int):
        run_tag = f'{run_id}'
    else:
        run_tag = f'{run_id[0]:0.0f}_{run_id[-1]:0.0f}'
        run_id = ','.join([f'{i:0.0f}' for i in run_id])

    # define timetrace path #
    trace_ascii = os.path.join(dirs['output'],
                               f'timetracei_{run_tag}.dat')
    trace_hdf5 = trace_ascii.replace('.dat', '.h5')
    with hf.File(trace_hdf5, 'w') as hf_:
        hf_.create_dataset('var indices', data=var_idx)
        hf_.create_dataset('kx indices', data=kx_idx)
        hf_.create_dataset('ky indices', data=ky_idx)

    # write command line parameters file #
    cl_dict = {'data path': dirs['problem'],
               'output path': dirs['output'],
               'start time': times[0],
               'end time': times[-1],
               'sparse factor': 1,
               'run ID': run_id,
               'species': 0,
               'diag': 'timetrace'}
    write_cl_diag(cl_dict)

    # read in timetrace diagnostic #
    diag_dir = get_diag_dir()
    tt_diag = os.path.join(diag_dir, 'prog', 'timetrace.pro')
    with open(tt_diag, 'r', encoding='utf-8') as file:
        diag_bakup = file.readlines()
    diag_new = diag_bakup.copy()

    # loop over timetrace variables #
    wdir = os.getcwd()
    os.chdir(diag_dir)
    for i, var in enumerate(var_idx):
        diag_new[31] = diag_bakup[31].replace('e.var = 0',
                                              f'e.var = {var}')
        for j, kx in enumerate(kx_idx):
            diag_new[32] = diag_bakup[32].replace('e.kxind = 0',
                                                  f'e.kxind = {kx:0.0f}')
            for k, ky in enumerate(ky_idx):
                nline = diag_bakup[33]
                diag_new[33] = nline.replace('e.kyind = par.ky0_ind EQ 0',
                                             f'e.kyind = {ky:0.0f}')

                # write temporary diagnostic file #
                with open(tt_diag, 'w') as file:
                    for line in diag_new:
                        file.write(line)

                # execute IDL diagnostic in child process #
                proc_idl = subprocess.Popen(shlex.split('idl cl_diag.pro'),
                                            stdin=subprocess.PIPE,
                                            stdout=subprocess.PIPE,
                                            stderr=subprocess.PIPE,
                                            text=True)
                proc_out, proc_err = proc_idl.communicate(input='exit')
                if proc_err:
                    raise RuntimeError(f"IDL call failed in timetrace (var={var}, kx={kx}, ky={ky}).\nSTDERR:\n{proc_err}")

                # read in timetrace data #
                with open(trace_ascii, 'r') as file:
                    lines = file.readlines()
                    for l, line in enumerate(lines):
                        if line[0] != '#':
                            ibeg = l
                            break
                    if i == 0 and j == 0 and k == 0:
                        nt = len(lines[ibeg::])
                        time = np.empty(nt)
                        fluc = np.empty(nt, dtype=complex)
                        for l, line in enumerate(lines[ibeg::]):
                            line_split = line.strip().split()
                            arr = np.array([float(x) for x in line_split])
                            time[l] = arr[0]
                            fluc[l] = arr[1]+1j*arr[2]
                    else:
                        fluc = np.empty(nt, dtype=complex)
                        for l, line in enumerate(lines[ibeg::]):
                            line_split = line.strip().split()
                            arr = np.array([float(x) for x in line_split[1::]])
                            fluc[l] = arr[0]+1j*arr[1]

                # write out to HDF5 file #
                with hf.File(trace_hdf5, 'a') as hf_:
                    hf_.create_dataset(f'var={var}, kx={kx}, ky={ky}', data=fluc)
                    if i == 0 and j == 0 and k == 0:
                        hf_.create_dataset('time', data=time)
    os.chdir(wdir)

    # reinstantiate original diagnostic file #
    with open(tt_diag, 'w') as file:
        for line in diag_bakup:
            file.write(line)
