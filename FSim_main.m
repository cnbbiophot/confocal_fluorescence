% Main program that launches the simulation of a confocal measurement of
% diffusive particles. It simulates a confocal microscope with one pixel
% with gaussian focus. The input variables is written in this main code and
% the input data is a struct varible with the following fields
%
% The struct has the name coodData[number_of_slice] and can be slices of
% the whole simulation (recommendation: no more than 1E5 steps and 1000 
% particles per slice, running from 1 to numbes_slices)
%
% Fields in struct:
% num_atoms [number of atoms in simulation]
% num_frames [number of time frames in the whole simulation]
% time_step [in ps]
% trajectory : double array with dimensions [num_particles, 3, num_steps_slice]

% The program is coded for a cube box

% Parameters (written inside the code):
% XMIN [um] origin of coordinales of the box
% XMAX [um]
% w0 [um]
% z0 [um]
% bin_size [s] // binning of the signal (has to be bigger than simulation binning)
% photon_mean [ph/s]
% photons_BK [ph/s] not used so far
% ispython % If the variable comes from python simulations
% several_species % If there are several species in the simulation

% Filenames : 
% name_f % String with the name of the gromacs variable folder
% name_simVar % String with the name of the gromacs variable
% name_pyf_folder % String with the name of the python variable folder
% name_pyf % String with the name of the python variable

% agv, 21oct2020
% agv, 27feb2021 ordered and set to publish
%% Parameter values
clear all

XMIN = 0; % in um
XMAX = 5;

w0 = 0.2; % in um %For calibration 200305
z0 = 1.3;
bin_size = 1e-5; % in s % binning (has to be bigger than simulation bin)

photon_mean = 1e4; % Number of photons per s (from experiment) % from experiment: 6672
photons_BK = 0; % Number of photons in background % For the moment we omit background

%% PARAMETER CHECK and OPTIONS

ispython = false; % When imputing data from python simulations (csv format)
several_species = false;
isbig = true; % If it has already sliced the signal in several coodDataX or python variables
Select_Distribution = true;

%% Filenames

name_f = 'D:\Users\Arturo\Matlab_sims\gromacs_variables\all_species\';
name_simVar = '21_400LUV_b5um_dt10us_t10s';

name_pyf = 'PY_p2000_b5000nm_D90_dt0.010000ms_t10000ms';
name_pyf_folder = 'D:\Users\Arturo\Python\trajectory_simulation\';

name_f_save = 'D:\Users\Arturo\Matlab_sims\fluorescence_simulations\210219_PCHlimits_1e4ph\';

%% Repeat experiment number_repetitions times

number_repetitions = 3;
for jk = 1:number_repetitions

do_experiment = true; % call several files to perform whole titration
[num_particles_exp, desired_number_part_exp, distribution_fluo_poiss,...
        k_average_exp, name_pyf_exp_dum, name_simVar_exp_dum, name_save_Var_exp_dum] = FSim_name_variables_and_sim_para();

if do_experiment
    i_exper_max = length(desired_number_part_exp); % number of points in titration
else % not do_experiment
    name_save_Var_exp_dum = false;
    i_exper_max = 1; % only one point in simulation
end

% Repeat the simulation for every point in titration
for i_exper = 1:i_exper_max

if do_experiment
    if ispython
        name_pyf = name_pyf_exp_dum{i_exper};
    else
        name_simVar = name_simVar_exp_dum{i_exper};
    end
    
    name_save_Var = name_save_Var_exp_dum{i_exper}; % name to save this titration point
    desired_number_part = desired_number_part_exp(i_exper); % desired number of particles randomly chosen this titration point
    num_particles = num_particles_exp(i_exper); % total number of particles in parent simulation input data in this titration point
    k_average = k_average_exp(i_exper); % average number of fluorophores per particle in this titration point
end
    
%% Simulation process
% Section to input number of fluorophores and prob distributions
if Select_Distribution; warning('number of particles, fluorophores and probability distribution chosen by the user. Make sure it fit with the variable parameters'); end

if Select_Distribution
    num_species_total = 15; % max number of allowed fluorophores
    num_fluoroph_homog = 1;
%     warning('I put 2 fluorophores per LUV');
    
    if do_experiment
        if distribution_fluo_poiss(i_exper)
            [ProbDist_part , Num_Fluorophores] = FSim_compute_distribution_fluorophores(num_particles, 'boltzmann', num_species_total, k_average, desired_number_part);
        else
            % I use the k_average as the number of fab in each distribution
            [ProbDist_part , Num_Fluorophores] = FSim_compute_distribution_fluorophores(num_particles, 'homogeneaus_number', desired_number_part, num_fluoroph_homog);
        end
    else % not do_experiment
        desired_number_part = 300;
        k_average = 2.07;
        num_particles = 400; % if only one file is loaded
        distribution_fluo_poiss = 0; % 1 if it is poisson
        
        if distribution_fluo_poiss
            [ProbDist_part , Num_Fluorophores] = FSim_compute_distribution_fluorophores(num_particles, 'poisson', num_species_total, k_average, desired_number_part);
        else
            [ProbDist_part , Num_Fluorophores] = FSim_compute_distribution_fluorophores(num_particles, 'homogeneaus_number', desired_number_part, num_fluoroph_homog);
        end
    end
    
else % Select_Distribution
    ProbDist_part = ones(num_particles,1); % Vertical vector with the number of fluorophores
    Num_Fluorophores = ones(num_particles,1);
end

    FSignal = FSim_set_data_for_simulation(XMIN, XMAX, w0, z0, bin_size, photon_mean, photons_BK,...
        ispython, several_species, isbig, name_f, name_simVar, name_save_Var, name_pyf_folder,...
        name_pyf, name_f_save, ProbDist_part, Select_Distribution, Num_Fluorophores);


end

end

disp('Done whole simulations set!')

