import os
import glob 
import shutil 
import argparse 
import GPUtil

import numpy as np
from mpi4py import MPI

import simtk.unit as u

from MD_utils.openmm_simulation import openmm_simulate_relaxation_exp # openmm_simulate_amber_npt 

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

print(f"{rank} of {size} printed") 

parser = argparse.ArgumentParser()
# parser.add_argument("-f", "--pdb_file", dest="f", help="pdb file")
# parser.add_argument("-p", "--topol", dest='p', help="topology file")
# parser.add_argument("-c", help="check point file to restart simulation")
parser.add_argument("-l", "--length", default=10, help="how long (ns) the system will be simulated")
# parser.add_argument("-g", "--gpu", default=0, help="id of gpu to use for the simulation")
args = parser.parse_args() 

deviceIDs = GPUtil.getAvailable(limit=size, maxLoad = 0.5, maxMemory = 0.5)
n_gpu = len(deviceIDs) 

input_paths = sorted(glob.glob('..//md_setup/complex_sim/examples/input_*_gmx'))
n_sim = len(input_paths)
# print(n_sim, input_paths) 
# exit()
print(n_sim, input_paths[rank]) 

proj_dir = os.getcwd() 

for i in range(rank, n_sim, size): 
    gpu_index = deviceIDs[rank]  # args.gpu # os.environ["CUDA_VISIBLE_DEVICES"]

    input_path = input_paths[i]
    input_path = os.path.abspath(input_path) 
    run_label = os.path.basename(input_path).replace('input_', '') 
    pdb_file = os.path.abspath(input_path + f'/{run_label[:-4]}.pdb')
    top_file = os.path.abspath(input_path + f'/{run_label[:-4]}.top')

    if os.path.exists(pdb_file) and os.path.exists(top_file): 
        print(rank, pdb_file, top_file) 
    else: 
        print("Error: ", input_path) 

    sim_path = f"run_{run_label}"
    os.makedirs(sim_path, exist_ok=True) 
    os.chdir(sim_path) 
    shutil.copy2(pdb_file, f"{run_label}.pdb") 
    shutil.copy2(top_file, f"{run_label}.top") 
    print(f"running simulation at {sim_path} on {gpu_index}") 

    openmm_simulate_relaxation_exp(pdb_file, top_file,
                             GPU_index=gpu_index,
                             temperature=310,
                             output_traj="output.dcd",
                             output_log="output.log",
                             # output_cm='output_cm.h5',
                             # ff_setup='charmm',
                             # relaxation=True, 
                             # constraint_group='UNK',
                             report_time=50*u.picoseconds,
                             sim_time=float(args.length)*u.nanoseconds)

    continue
    # check_point = None
    try: 
        openmm_simulate_relaxation_exp(pdb_file, top_file,
                                 GPU_index=gpu_index,
                                 temperature=310,
                                 output_traj="output.dcd",
                                 output_log="output.log",
                                 # output_cm='output_cm.h5',
                                 # ff_setup='charmm',
                                 # relaxation=True, 
                                 # constraint_group='UNK',
                                 report_time=50*u.picoseconds,
                                 sim_time=float(args.length)*u.nanoseconds)
    except Exception as e:  
        print(f"simulation on {sim_path} failed", e) 

    print(f"simulation done at gpu {gpu_index}") 
    os.chdir(proj_dir) 
