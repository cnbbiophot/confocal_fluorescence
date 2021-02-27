function final_FSignal = FSim_fluorescence_simulation_from_traj...
    (XMIN, XMAX, w0, z0, photon_mean, photons_BK, deltaT_SIM, bin_size_bin,...
    coodData, several_species, ProbDist_part, Select_Distribution, Num_Fluorophores)

% Simulation of confocal fluorescence (gaussian focus) from trajectory
% simulation variable computed with gromacs

% The concept follows, particle by particle:
% 0- Takes a whole time simulation from a single particle
% 1- Check which parts of the trajectory are in focus 
% 2- Compute value of gaussian focus at those points
% 3- Compute photon emmision with/without Poisson noise
% 4- Add the photons to the bin of the general signal

% coodData is a strut ->
% (e.g.) num_atoms: 20, num_frames: 500001, time_step: 5000000, trajectory: [20×3×500001 double]

% deltaT_GROM is the size of the temporal step of the gromacs simulation. All the particles move at
% once. E.g. 5us
% bin_size_bin is the size of the time step of the desired binning

% agv, 21sept2020
% agv, 6oct2020 new concept of fluorescence simulation 
% agv, 21oct2020 fix problems

dum_coodData_all = struct2cell(coodData);
num_slices = size(dum_coodData_all,3); % See the number of slices in the signal

p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)
   parpool(4)
end

for i_slices = 1:num_slices % go slice by slice
       
    if num_slices > 1
        disp(['Compute slice ' num2str(i_slices) ' out of ' num2str(num_slices)])
    end
    
    % Split the data into different variables
    dum_coodData = dum_coodData_all(:,:,i_slices);
    num_particles = cell2mat(dum_coodData(1));
    num_particles_total = sum(num_particles);
    
    % Generate partProp, a vector with information of each particle 
%     partProp(:,1) = 1:num_particles_total; % name of the particles
%     partProp(:, 2) = Num_Fluorophores; % number of fluorophores of each particle
%     partProp(:, 3) = ProbDist_part; % set whether the particle is silent or contribute to the signal
    
    if several_species % If there are several species, the variables are differently distributed
        sim_steps = round(cell2mat(dum_coodData(3)) / num_slices); % Number of total simulation step
        num_species = length(num_particles);
        num_fluoroph = cell2mat(dum_coodData(2));
        for i = 1:num_species
            dum_num_P = cumsum(num_particles);
            if i == 1
                s = 1;
            else
                s = 1 + dum_num_P(i-1);
            end
            e = dum_num_P(i);
            % The particles are ordered as this [S1, S1, S1, S2, S2, S2]
            % and then they are grupped
            partProp(s:e,2) = ones(num_particles(i),1) * num_fluoroph(i);
        end
        PosParticlest = cell2mat(dum_coodData(5)); % with several species the trajectories are in 5th position
    else % single species
        partProp = zeros(num_particles_total, 2); % Properties of the particles [numberP #fluorophores]
        partProp(:,2) = ones(num_particles_total,1); % Se the number of fluorophores per LUV or Free fab
        % In this case, only one fluorophore per particle (only one species)
        partProp(:,1) = 1:num_particles_total; % name of the particles
        PosParticlest = cell2mat(dum_coodData(4));
        % % sim_steps = round(cell2mat(dum_coodData(2)) / num_slices); % Number of total simulation step
        sim_steps = length(PosParticlest(1,1,:)); % If the data comes sliced but not with all the pieces at once (gromacs + isbig)
    end

    max_time_sim = deltaT_SIM * sim_steps;
    num_bins = round(max_time_sim / bin_size_bin);
    FSignal = zeros(num_bins,2); % Binned signal to input photon counts
    FSignal(:,1) = double(((1:num_bins) - 1)).*bin_size_bin;
    
    % To gather all the slices, create the main trace only at the beggining
    if ~exist('final_FSignal')
        final_FSignal = zeros(num_bins*num_slices,2); % Binned signal to input photon counts
        final_FSignal(:,1) = double(((1:num_bins*num_slices) - 1)).*bin_size_bin;
    end

%     f = waitbar(0,'Compute photon trace');
    
    tic
    deltaT_BK = bin_size_bin; % set tipical time interval
    centre_PSF = (XMAX - XMIN)/2;
    
    if Select_Distribution == true
        partProp(:, 3) = ProbDist_part;
        partProp(:, 2) = Num_Fluorophores;
    else
        partProp(:, 3) = ones (length(partProp(:, 2)),1);
    end
    
    PosParticlest = PosParticlest * 1e-3; % from nm to um (because of GROMACS data
    
    % Compute photons particle by particle
    for ipart = 1:num_particles_total 
        
        if partProp(ipart, 3) == 0; continue; end % if the particle is silent, don't contribute to signal
        
        % re-writte some variable for parfor loop (we have to overwrite totally)
        
        % Find the positions that lie in the focus
        focus_pos = isinfocus(PosParticlest(ipart, :, :), w0, z0, centre_PSF, sim_steps); % see when they are in focus
%         F = parfeval(@isinfocus, 1, PosParticlest(ipart, :, :), w0, z0, centre_PSF, sim_steps); % this is slower
%         focus_pos = fetchOutputs(F);
        
        focus_pos = squeeze(PosParticlest(ipart, :, :)) .* [focus_pos; focus_pos; focus_pos]; % leave only nonzero those that are in focus
        
        FocusValues = exp(-2*((focus_pos(1,:)-centre_PSF).^2 + (focus_pos(2,:)-centre_PSF).^2) ./ w0^2 ...
            - 2*(focus_pos(3,:)-centre_PSF).^2 ./ z0^2); % Values of the focus at particle positions

        % Compute poisson emission
        pois_mean_inSIMBIN = FocusValues * photon_mean * deltaT_SIM * partProp(ipart, 2);
        photons_em_step = poissrnd(pois_mean_inSIMBIN); % With Poisson randomness THIS IS NEEDED TO HAVE A REAL SIGNAL
        
        binnedFS = zeros(1,num_bins); % reset counting of photons for every particle
        % binnedFS has to be there in order to use parfor
        parfor ibin = 1:num_bins % Binning of the photon signal
            photons_inbin = photons_em_step(1 + round((ibin - 1) * sim_steps / num_bins) ...
                : round((ibin) * sim_steps / num_bins));
            binnedFS(ibin) = binnedFS(ibin) + sum(photons_inbin);
        end

        FSignal(:,2) =  FSignal(:,2) + binnedFS' ;
        
%     waitbar(ipart/num_particles_total,f)            
    end 
    
    if photons_BK ~= 0 % if there is noise in the system
        add_noise = poissrnd(photons_BK*deltaT_SIM, num_bins, 1); % Generate random Poisson noise through the signal
        FSignal(:,2) = FSignal(:,2) + add_noise;
    end
    
    toc
%     close(f)
   
    final_FSignal(1 + (i_slices-1) * num_bins : i_slices * num_bins,2) = FSignal(:,2);
    
end

% delete(p) % Stop parpool    
end

function focus_time_vector = isinfocus(positions, w0, z0, centre_PSF, sim_steps)
% Function that looks into the positions and see which of them are within a
% box of size [2*w0, 2*w0, 2*z0]. Return a vector with 1 (in focus) and 0 for every
% time of the trajectory

LIMupX = centre_PSF + w0;
LIMdownX = centre_PSF - w0;
LIMupZ = centre_PSF + z0;
LIMdownZ = centre_PSF - z0;

positions = squeeze(positions);

outX = find (positions(1,:) > LIMdownX & positions(1,:) < LIMupX);
outY = find (positions(2,:) > LIMdownX & positions(2,:) < LIMupX);
outZ = find (positions(3,:) > LIMdownZ & positions(3,:) < LIMupZ);

dum_infocus = zeros(1,sim_steps); % Sparse does not work here

infocusX = dum_infocus;
infocusX(outX) = 1;
infocusY = dum_infocus;
infocusY(outY) = 1;
infocusZ = dum_infocus;
infocusZ(outZ) = 1;

focus_time_vector = infocusX.*infocusY.*infocusZ; % keep only positions where the three coordinates are in focus

end
