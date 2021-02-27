function [ProbDist_part , Num_Fluorophores] = FSim_compute_distribution_fluorophores(num_particles, ...
    type_distribution, varargin)

% Computes the distribution of fluorophores for the given condition
%
% ProbDist_part: a (num_particles, 1) vector woth ones and zeros that
%       allows to avoid the computation of one or more particles as if they
%       where silent (dark)
% Num_Fluorophores: a (num_particles, 1) vector that contains the number of
%       fluorophores of each particle
%
% agv, 04Nov2020
% agv, 12feb2021 insert Boltzmann distribution

    if strcmp(type_distribution, 'poisson')
        
        num_species_max = varargin{1}; % max number of allowed fluorophores
        k_average = varargin{2};
        
        if nargin == 5
            num_particles_dum = varargin{3};  
        else
            num_particles_dum = num_particles;
        end
        
        dist = poisspdf(0:num_species_max, k_average);
        dist = floor(dist*num_particles_dum);
        
        for i_sp = 1:num_species_max
            indexes = cumsum(dist);
            if i_sp == 1
                Num_Fluorophores(1:indexes(i_sp)) = 0;
            elseif dist(i_sp) ~= 0
                Num_Fluorophores(indexes(i_sp - 1) + 1 : indexes(i_sp)) = i_sp - 1;  
            else
                continue;
            end
        end
        
        if length(Num_Fluorophores) < num_particles
            Num_Fluorophores(length(Num_Fluorophores)+1 : num_particles) = 0;
        end
        
        Num_Fluorophores = Num_Fluorophores';
        
        indices = randperm(num_particles, num_particles);
        Num_Fluorophores(indices) = Num_Fluorophores; % randomize the particles chosen
        
        ProbDist_part = ones(num_particles, 1);
        ProbDist_part(Num_Fluorophores==0) = 0;
        
    elseif strcmp(type_distribution, 'boltzmann')
        
        num_species_max = varargin{1}; % max number of allowed fluorophores
        k_average = varargin{2};
        
        lambda = k_average; % Boltzmann distribution is defined with Lambda
        
        if nargin == 5
            num_particles_dum = varargin{3}; % number of particles in the variable to simulate
        else
            num_particles_dum = num_particles;
        end
        
        dist = zeros(1,num_species_max+1);
        n = 0:num_species_max;
        dist = (lambda .^ n) ./ ( (lambda + 1).^(n+1) );
        dist = round(dist*num_particles_dum);
        
        for i_sp = 1:num_species_max % Assign each particle the number of fluorophores
            indexes = cumsum(dist);
            if i_sp == 1
                Num_Fluorophores(1:indexes(i_sp)) = 0;
            elseif dist(i_sp) ~= 0
                Num_Fluorophores(indexes(i_sp - 1) + 1 : indexes(i_sp)) = i_sp - 1;  
            else
                continue;
            end
        end
        
        if length(Num_Fluorophores) < num_particles
            Num_Fluorophores(length(Num_Fluorophores)+1 : num_particles) = 0;
        end
        
        Num_Fluorophores = Num_Fluorophores';
        
        indices = randperm(num_particles, num_particles);
        Num_Fluorophores(indices) = Num_Fluorophores; % randomize the particles chosen so every distribution is different
        
        ProbDist_part = ones(num_particles, 1);
        ProbDist_part(Num_Fluorophores==0) = 0; % make 0 all the rest of particles so they don't undergo fluorescence
        
        
    elseif strcmp(type_distribution, 'homogeneaus_number')
        
        desired_number_part = varargin{1};
        number_fluoroph_homo = varargin{2};
        
        if desired_number_part >= num_particles
            disp('The desired number of particles has to be smaller than ')
        end
        
        Num_Fluorophores = ones(num_particles, 1)*number_fluoroph_homo;
        
        Num_Fluorophores(1 : num_particles-desired_number_part) = 0;
        indices = randperm(num_particles, num_particles);
        Num_Fluorophores(indices) = Num_Fluorophores; % randomize the particles chosen
        
        ProbDist_part = Num_Fluorophores;    
        
    elseif strcmp(type_distribution, 'exponential')
        
        num_species_max = varargin{1}; % max number of allowed fluorophores
        k_average = varargin{2};
        
        if nargin == 5
            num_particles_dum = varargin{3};  
        else
            num_particles_dum = num_particles;
        end
        
        dist = exppdf(0:num_species_max, k_average);
        dist = dist/sum(dist); % NORMALIZE
        dist = floor(dist*num_particles_dum);
        
        for i_sp = 1:num_species_max
            indexes = cumsum(dist);
            if i_sp == 1
                Num_Fluorophores(1:indexes(i_sp)) = 0;
            elseif dist(i_sp) ~= 0
                Num_Fluorophores(indexes(i_sp - 1) + 1 : indexes(i_sp)) = i_sp - 1;  
            else
                continue;
            end
        end
        
        if length(Num_Fluorophores) < num_particles
            Num_Fluorophores(length(Num_Fluorophores)+1 : num_particles) = 0;
        end
        
        Num_Fluorophores = Num_Fluorophores';
        
        indices = randperm(num_particles, num_particles);
        Num_Fluorophores(indices) = Num_Fluorophores; % randomize the particles chosen
        
        ProbDist_part = ones(num_particles, 1);
        ProbDist_part(Num_Fluorophores==0) = 0;
        
        
    end