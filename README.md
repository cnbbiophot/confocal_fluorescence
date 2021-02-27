# confocal_fluorescence
This program simulates the fluorescence emmision of diffusing particles crossing a confocal focus. So far, the focus has a Gaussian shape.

## FSim_main.m

The main program  is not a matlab function but a bunch of code that contains the parameters and the initial calling to the simulation function. The input variables is written in this main code and the input data is a struct varible with the following fields:

Fields in struct (follow the same order):
1. num_atoms [number of atoms in simulation]
2. num_frames [number of time frames in the whole simulation]
3. time_step [in ps]
4. trajectory : double array with dimensions [num_particles, 3, num_steps_slice]

This struct has the name coodData[number_of_slice] and can be sliced. My recommendation to avoid memory problems is that slices have no more than 1000 particles and 1E5 steps.

The program is coded for data diffusing in a cubic box.

Some parameters to choose are: 
* Photon counts per molecule per second
* Size of focus to analyse
* Photon background counts
* Type of fluorophore distributions if loaded particles are used

## FSim_compute_distribution_fluorophores.m

Configuration function of names and simulation parameters to use in FSim

Parameters used are:
1. num_particles_exp is the number of particles in each condition
2. desired_number_part_exp is the number of particles that will have non zero number of fluorophores
3. distribution_fluo_poiss logical variable that tell whether ther is an specific fluorophore distribution
4. k_average_exp average number of fluorophores pwe condition, 0 if there is no explicit distribution
5. name_simVar and name_pyf_dum are cell variables with the gromacs/python variables names

## FSim_fluorescence_simulation_from_traj.m

function [ProbDist_part , Num_Fluorophores] = FSim_compute_distribution_fluorophores(num_particles, ...
    type_distribution, varargin)

Computes the distribution of fluorophores for the given condition in every particle. When less than all the available number of particles are used, the function chooses randomli among the particles to generate as much variability as possible in each simulation.

Output variables:  
1. ProbDist_part: a (num_particles, 1) vector with ones and zeros that allows to avoid the computation of fluorescence of one or more particles if they are silent (dark)
2. Num_Fluorophores: a (num_particles, 1) vector that contains the number of fluorophores of each particle

## FSim_name_variables_and_sim_para.m

Configuration variable with the information of the variables to input and the desired parameters, and the automatic generator of output names depending on the pourpose.

## FSim_set_data_for_simulation.m

Function that put the variables (either gromacs of python format) in the correct shape and format to perform the simulation.

## compute_save_brownian_trajectory.py 

Simple python program that computes and save simple brownian motion in 3D and stores the trajectoryes in csv format.
