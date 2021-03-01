# Confocal microscopy fluorescence simulation

This program simulates the fluorescence emmision of diffusing particles crossing a confocal focus. So far, the focus has a Gaussian shape. The program accepts background noise and assumes a poisson emission-detection model for photons.

# The fluorescence simulation program

## FSim_main.m

The main program  is a matlab script that contains the main parameters and the initial calling to the simulation function. 

The program is coded for simulated traces of particles diffusing in a cubic box, but it can be easily rewritten to have other shapes of simulation boxes.

Some parameters to choose are: 
* Photon counts per molecule per second
* Size of focus to analyse
* Background photon counts
* Type of fluorophore distributions if loaded particles are used

## FSim_compute_distribution_fluorophores.m

Configuration function of names and simulation parameters to use in FSim.

Parameters used are:
1. num_particles_exp is the number of particles in each condition
2. desired_number_part_exp is the number of particles that will have non zero number of fluorophores
3. distribution_fluo_poiss logical variable that tell whether ther is an specific fluorophore distribution
4. k_average_exp average number of fluorophores pwe condition, 0 if there is no explicit distribution
5. name_simVar and name_pyf_dum are cell variables with the gromacs/python variables names

## FSim_fluorescence_simulation_from_traj.m

Computes the distribution of fluorophores for the given condition in every particle. When less than all the available number of particles are used, the function chooses randomly among the particles to generate as much variability as possible in each simulation.

## FSim_name_variables_and_sim_para.m

Configuration variable with the information of the variables to input and the desired parameters, and the automatic generator of output names depending on the pourpose.

## FSim_set_data_for_simulation.m

Function that put the variables (either gromacs of python format) in the correct shape and format to perform the simulation. Also handles the different slicing of the signal to simulate slice by slice and then stick toguether all the signals.

## FSim_fluorescence_simulation_from_traj.m

Main simulation program. It assumes a poisson emmission-detection process for each fluorophore. It extracts the parts of the trajectory that are in focus with a sufficient tolerant limits. Then computes the photon counts for each fluorophore separately and add each fluorophore emmision to the final signal with the correct binning.

# Input data and diffusive simulations

Normally, we use input data from a gromacs diffusion simulation. The input data is a matlab struct varible with the following fields:

1. num_atoms [number of atoms in simulation]
2. num_frames [number of time frames in the whole simulation]
3. time_step [in ps]
4. trajectory : double array with dimensions [num_particles, 3, num_steps_slice]

This struct has the name *coodData[number_of_slice]* (e.g. coodData1, coodData23, etc.) and can be sliced. My recommendation to avoid memory problems is that slices have no more than 1000 particles and 1E5 steps. Otherwise, the variables are too big to be handled by a 16Gb RAM memory and Matlab can crash. You can use as many slices as you need.

## compute_save_brownian_trajectory.py 

It is a simple python program that computes and save simple brownian motion in 3D and stores the trajectoryes in csv format. When used with FSim_main.m remember to assign 
```
ispython = true
```
