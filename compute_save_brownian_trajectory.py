# Dependencies
import numpy
import matplotlib
from matplotlib import pylab, mlab, pyplot
np = numpy
plt = pyplot

from IPython.display import display
from IPython.core.pylabtools import figsize, getfigs

from pylab import *
from numpy import *

import matplotlib.pyplot as plt
import time
from scipy.stats import norm
import copy
import csv

# agv, 20oct2020


def brownian_only_numpy(total_sim_time, time_step, num_of_mol, D, width, height, depth, track_arr_init):
    """

    Inputs:
    total_simulation_time: Total simulation time in ms.
    time_step:             The duration of each time step ms.
    num_mol:               The number of molecules in the simulation.
    D:                     The diffusion rate.
    width:                 The width of the simulation.
    height:                The height of the simulation.
    Outputs:
    track_arr:             A dictionary where each track number (e.g. track_arr[0]) contains the track data [0,:] [1,:]
    """

    # Number of steps.
    num_of_steps = int(round(float(total_sim_time) / float(time_step), 0))

    print('num_of_steps', num_of_steps)
    # Calculates length scales
    scale_in = np.sqrt(2.0 * (float(D) * 1e3) * float(time_step))

    if track_arr_init == 0:
        # Randomly generates start locations
        start_coord_x = (np.random.uniform(0.0, 1.0, num_of_mol)) * width
        start_coord_y = (np.random.uniform(0.0, 1.0, num_of_mol)) * height
        start_coord_z = (np.random.uniform(0.0, 1.0, num_of_mol)) * depth
    else:
        # select initial positions
        start_coord_x = np.zeros((num_of_mol,1))
        start_coord_y = np.zeros((num_of_mol,1))
        start_coord_z = np.zeros((num_of_mol,1))
        # set coordinates as last and move them so they don't repeat points
        for particle in range(num_of_mol):
            start_coord_x[particle] = track_arr_init[particle][0] + norm.rvs(size=[1]) * scale_in
            start_coord_y[particle] = track_arr_init[particle][1] + norm.rvs(size=[1]) * scale_in
            start_coord_z[particle] = track_arr_init[particle][2] + norm.rvs(size=[1]) * scale_in

    track_arr = {}

    # This can be done as one big matrix, but can crash system if large so I break it up into the individual molecules.
    for b in range(0, num_of_mol):
        print('processing tracks: ', np.round((float(b) / float(num_of_mol)) * 100, 1), '%')
        track_arr[b] = np.zeros((3, num_of_steps))
        track_arr[b][0, 0] = start_coord_x[b]
        track_arr[b][1, 0] = start_coord_y[b]
        track_arr[b][2, 0] = start_coord_z[b]
        rand_in = norm.rvs(size=[3, num_of_steps]) * scale_in
        track_arr[b][:, 1:] += rand_in[:, 1:]
        track_arr[b] = np.cumsum(track_arr[b], 1)
        out = track_arr[b]
        mod = np.zeros((out.shape))
        mod[0, :] = np.floor(track_arr[b][0, :].astype(np.float64) / height)
        mod[1, :] = np.floor(track_arr[b][1, :].astype(np.float64) / width)
        mod[2, :] = np.floor(track_arr[b][2, :].astype(np.float64) / depth)
        track_arr[b] = np.array(out - ([mod[0, :] * height, mod[1, :] * width, mod[2, :] * depth]))

    return track_arr

def compute_and_store_trajectory(width,height,depth,total_sim_time,num_of_mol,time_step,D,max_number_steps_piece):

    num_of_steps = int(round(float(total_sim_time) / float(time_step), 0))

    # trajectories (separate variable)
    name_folder = r"D:\\Users\\Arturo\\Python\\SimulatedTrajectories\\"
    name_file_traj = r"PY_p%d_b%dnm_D%d_dt%fms_t%dms.csv" % (num_of_mol, width, D, time_step, total_sim_time)
    name_file_param = r"SET_" + name_file_traj

    name_file_traj = name_folder + name_file_traj
    name_file_param = name_folder + name_file_param

    with open(name_file_param, 'w') as f:
        f.write("%f\n" % (num_of_mol))
        f.write("%f\n" % (num_of_steps))
        f.write("%f\n" % (time_step))

    t_total = time.time()
    #Simulates the brownian motion.

    if num_of_steps > max_number_steps_piece:
        number_of_pieces = round(num_of_steps / max_number_steps_piece)
        print('Computing trajectory in ', number_of_pieces, ' pieces.')
        num_of_steps_piece = max_number_steps_piece
        total_sim_time = round(total_sim_time / number_of_pieces)
    else:
        number_of_pieces = 1
        num_of_steps_piece = num_of_steps

    for piece in range(number_of_pieces):

        t_piece = time.time()

        # compute trajectories
        if piece == 0:
            track_arr = {}
            track_arr[piece] = brownian_only_numpy(total_sim_time, time_step, num_of_mol, D, width, height, depth, 0)
        else:
            track_arr_init_dum = {}
            for particle in range(num_of_mol):
                track_arr_init_dum[particle] = np.zeros((3, 1))
                dum1 = track_arr[piece-1]
                dum2 = dum1[particle]
                dum3 = dum2[:, -1] # get the last elememt
                track_arr_init_dum[particle] = dum3

            del track_arr, track_arr_dum

            # track_arr is a dictionary and we refer to each piede by its name (0, 1, 2) no position
            track_arr = {}
            track_arr[piece] = brownian_only_numpy(total_sim_time, time_step, num_of_mol, D, width, height, depth, track_arr_init_dum)
        # Save trajectories
        if number_of_pieces > 1:
            name_file_traj_dum = name_file_traj[:-4] + "_piece%d" % (piece + 1) + ".csv"
            print(name_file_traj_dum)
        else:
            name_file_traj_dum = name_file_traj

        # Prepare files to save
        track_arr_dum = {}
        # to treat unhashable arrays
        for particle in range(num_of_mol):
            track_arr_dum[particle] = np.zeros((num_of_steps_piece, 3))
            dum1 = track_arr[piece]
            dum2 = dum1[particle]
            dum3 = dum2[:,:]
            track_arr_dum[particle] = np.transpose(dum3)
            del dum1, dum2, dum3

        '''
        # Plot trajectory to check
        plt.plot(track_arr_dum[0][:10000,0], track_arr_dum[0][:10000,1])
        plt.axis([0, width, 0, depth])
        plt.show()
        '''

        # Store trajectories
        with open(name_file_traj_dum, 'w') as f:
            for particle in range(num_of_mol):
                print('Storing as CSV: ', np.round((float(particle) / float(num_of_mol)) * 100, 1), '%')
                for time_step_index in range(num_of_steps_piece):
                    f.write("%f,%f,%f\n" % (track_arr_dum[particle][time_step_index][0], track_arr_dum[particle][time_step_index][1], track_arr_dum[particle][time_step_index][2]))

        elapsed_P = time.time() - t_piece
        print('Elapsed time: ', elapsed_P)

    elapsed_T = time.time() - t_total
    print ('Elapsed time: ',elapsed_T)

    return None

width = 5000.0  # The width of the simulation. # Units in nm
height = 5000.0  # The height of the simulation. Should be greater than the active area radius. # Units in nm
depth = 5000.0  # Units in nm
total_sim_time = 10000  # Total time. ms (former value 30000)
num_of_mol = [2000]  # Number of molecules in simulation.
time_step = 0.01  # Time step in ms. (former value 0.04)
D = 90  # diffusion coefficient [um^2/s](former value 1.0)
max_number_steps_piece = 100000

number_sim = num_of_mol.__len__()

for i in range(number_sim):
    print('Simulation with ', num_of_mol[i], ' molecules')
    compute_and_store_trajectory(width, height, depth, total_sim_time, num_of_mol[i], time_step, D, max_number_steps_piece)
