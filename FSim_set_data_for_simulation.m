function FSignal = FSim_set_data_for_simulation(XMIN, XMAX, w0, z0, bin_size, photon_mean, photons_BK,...
    ispython, several_species, isbig, name_f, name_simVar, name_save_Var, name_pyf_folder,...
    name_pyf, name_f_save, ProbDist_part, Select_Distribution, Num_Fluorophores)

% Set the input variables with particle trajectories in the rigth shape to
% perform confocal microscopy fluorescence simulation

% XMIN [um] origin of coordinales of the box
% XMAX [um]
% w0 [um]
% z0 [um]
% bin_size s // binning of the signal (has to be bigger than simulation binning)
% photon_mean in ph/s
% photons_BK in ph/s not implemented yet
% ispython % If the variable comes from python simulations
% several_species % If there are several species
% name_f % String with the name of the gromacs variable folder
% name_simVar % String with the name of the gromacs variable
% name_pyf_folder % String with the name of the python variable folder
% name_pyf % String with the name of the python variable

% Python varibles will be in csv format and gromacs variables in matlab 
% struct format

% agv, 21oct2020
% agv, 27feb2021 get rid of simulations for matlab and RDME input data format

    disp(['Calculate traze for data with gromacs format'])
    if ispython 
        disp(['Calculated with python brownian simulation']) 
    else
        disp(['Calculated with gromacs'])
    end

    if several_species
        disp('***For SEVERAL species***')
    else
        disp('***For ONE species***')
    end
    name_simVar = {name_simVar};

    if ispython; name_simVar = {name_pyf}; end % set the name to the python name if it is the case

    DateString = datestr(now, 'yymmdd_HHMM');
    name_timeBIN = {[DateString '_' cell2mat(name_simVar) '_binned_']};

    name_f = {name_f}; % Transform the string into a cell with a string inside

    for i_names_loop = 1: length(name_simVar) % Display variables name

        name_var = {name_f name_simVar{i_names_loop}};
        name_var=strcat(name_var{1},name_var{2},name_var{3});
        disp(name_var{1})

    end
    
    if ispython

        tic
        disp('Loading simulation from python simulation in csv format')
        disp('Extracting data...')
%             struct_dum = csvread('C:\Users\singlemol\Documents\Python\trac_py.csv');

        % Open python setting variable
        name_pyf_setting = ['SET_' name_pyf];
        fid = fopen([name_pyf_folder name_pyf_setting '.csv'], 'rt');
        warning('Set properly the number of columns read if several species are input');
        set_struct_dum = cell2mat(textscan(fid,'%f%f%f','delimiter',','));
        fclose(fid);

        files_dum = dir([name_pyf_folder name_pyf '_piece*']); % check if the variable is in pieces
        number_pieces_traj = size(files_dum,1);

        if number_pieces_traj > 0 % Check 
            isbig = 1;
            disp(['Analyzing ' num2str(number_pieces_traj) ' pieces'])
        else
            isbig = 0;
            number_pieces_traj = 1;
        end

        for i_piece_t = 1:number_pieces_traj % Put the python data in gromacs format
            % Open trajectories
            if isbig % If isbig there will be many slices
                disp(['Analyze signal piece ' num2str(i_piece_t)])
                fid = fopen([name_pyf_folder files_dum(i_piece_t).name], 'rt');
            else
                fid = fopen([name_pyf_folder name_pyf '.csv'], 'rt');
            end
            struct_dum = cell2mat(textscan(fid,'%d%d%d','delimiter',',')); % read as int because is faster
            fclose(fid);

            disp('Extracted!')
            toc

            % build the struct as gromacs data            
            if several_species
                coodData.num_atoms = double(set_struct_dum(1,:)); % There will be several species
                coodData.num_fluoroph = double(set_struct_dum(2,:));
                coodData.num_frames = double(set_struct_dum(4,1));
                coodData.time_step = double(set_struct_dum(5,1)); % in ms                   
            else
                coodData.num_atoms = double(set_struct_dum(1,1));
                if isbig
                    coodData.num_frames = round(double(set_struct_dum(2,1)) / number_pieces_traj);
                else
                    coodData.num_frames = double(set_struct_dum(2,1));
                end
                coodData.time_step = double(set_struct_dum(3,1)); % in ms
            end

            simulation_bin = coodData.time_step * 1e-3; % is in ms

            % The trajectories are in order [traj_part1; traj_part2; etc]     

            num_frames = coodData.num_frames;
            num_atoms = sum(coodData.num_atoms); % Changed in 201021 TEST

            tic
            disp('Building variable')
            
            % Slice variable if it is too big (due to memory issues)
            frames_slice = 5e5; % proper 5e5

            if coodData.num_frames > frames_slice % perform slicing of the signal in pieces of 5s (only multiples of 5 s allowed)
                if mod(coodData.num_frames,frames_slice) ~= 0 ; error('The slicing of the signal has to be with even pieces'); end
                num_slices = coodData.num_frames/frames_slice; % For now, signals with only multiples of 5E5, it can be changed
                disp(['Slicing the signal in ' num2str(num_slices) ' pieces of ' num2str(frames_slice) ' steps'])
            else
                num_slices = 1;
                frames_slice = coodData.num_frames;
                coodData.trajectory = zeros(num_atoms, 3, coodData.num_frames);
            end

            % generate all the structures needed 
            for i_slices = 1:num_slices
                dum_coodData(i_slices)=coodData;
            end

            % Set the trajectory varibles in the right shape
            for i_slices = 1:num_slices
                if num_slices > 1 % this is very slow, but I don't know how to overcome it
                    disp(['Slice ' num2str(i_slices)]);
                    for i = 1:num_atoms
                    s = 1 + (i_slices-1)*frames_slice + (i-1) * num_frames;
                    e = s - 1 + frames_slice;
                    dum_coodData(i_slices).trajectory(i,:,:) = transpose(single(struct_dum(s:e,:)));
                    end
                else
                    for i = 1:num_atoms
                    s = 1 + (i_slices-1)*frames_slice + (i-1) * num_frames;
                    e = s - 1 + frames_slice;
                    coodData.trajectory(i,:,:) = transpose(single(struct_dum(s:e,:)));
                    end
                end                    
            end

             % gather all the variables in coodData, generate string for
             % all the pieces and evaluate it
             if num_slices > 1
                coodData = dum_coodData(1);
                for i_slices = 2:num_slices         
                    coodData = horzcat(coodData, dum_coodData(i_slices));
                end
             end

            toc
            disp('Finished building!')

            if isbig
                FSignal_dum{i_piece_t} = FSim_fluorescence_simulation_from_traj...
                            (XMIN, XMAX, w0, z0, photon_mean, photons_BK, simulation_bin, bin_size,...
                            coodData, several_species, ProbDist_part, Select_Distribution, Num_Fluorophores);
            else
                FSignal = FSim_fluorescence_simulation_from_traj...
                                (XMIN, XMAX, w0, z0, photon_mean, photons_BK, simulation_bin, bin_size,...
                                coodData, several_species, ProbDist_part, Select_Distribution, Num_Fluorophores);
            end
        end

        if isbig % gather all pieces of the signal
            FSignal = FSignal_dum{1};
            for i = 2:number_pieces_traj
                FSignal = vertcat(FSignal, FSignal_dum{i});
            end
            FSignal(:,1) = (1:length(FSignal(:,1))) * bin_size;
        end

    else % data format from gromacs simulation
        disp('Loading simulation from gromacs')
        disp('Extracting data...')
        
        if isbig % if the variable is too big to open one by one
            disp('As a big variable, we compute it piece by piece')
            variableInfo = who('-file', name_var{1}); % Check variables inside the variable loaded

            num_pieces = sum(strncmp(variableInfo,'coodData',8) == 1);
            disp(['Sliced in ' num2str(num_pieces) ' pieces']);

            for num_slice = 1:num_pieces
                disp(['Load and compute slice ' num2str(num_slice)])
                eval(['load(name_var{1}, ''coodData' num2str(num_slice) ''')']);

                eval(['coodData = coodData' num2str(num_slice) ';']); % name it coodData

                simulation_bin = coodData.time_step * 1E-12; % the value is in ps
                if bin_size < simulation_bin; error('The binning cannot be bigger than the simulation dT'); end

                % Perform fluorescence simulation
                FSignal_dum{num_slice} = FSim_fluorescence_simulation_from_traj...
                        (XMIN, XMAX, w0, z0, photon_mean, photons_BK, simulation_bin, bin_size,...
                        coodData, several_species, ProbDist_part, Select_Distribution, Num_Fluorophores);

                eval(['clear ''coodData' num2str(num_slice) '''']);
            end

            num_slices = num_slice;
            FSignal = vertcat(FSignal_dum{1}, FSignal_dum{2});

            if num_slices > 2
                for i = 3:num_slices
                    FSignal = vertcat(FSignal, FSignal_dum{i});
                end
            end

            FSignal(:,1) = (1:length(FSignal(:,1))) * bin_size;

        else % If the variable can be handled in the memory
            load(name_var{1})
            disp('Extracted!')
            if exist('coodData2') % Check the slicing and put the variables toguether
                warning(sprintf(['The input data is sliced.']))
                s=whos('coodData*');
                num_slices = size(s,1);
                coodData = coodData1;
                for i_slices = 2:num_slices         
                    eval(['coodData = horzcat(coodData, coodData' num2str(i_slices) ' );']);
                end
                simulation_bin = coodData1.time_step * 1E-12; % the value is in ps  
            else
                simulation_bin = coodData.time_step * 1E-12; % the value is in ps    
            end
                if bin_size < simulation_bin; error('The binning cannot be bigger than the simulation dT'); end

                FSignal = FSim_fluorescence_simulation_from_traj...
                        (XMIN, XMAX, w0, z0, photon_mean, photons_BK, simulation_bin, bin_size,...
                        coodData, several_species, ProbDist_part, Select_Distribution, Num_Fluorophores);
        end
    end


    if name_save_Var == false
        name_var = {name_f_save name_timeBIN {num2str(bin_size)} {'s.mat'}};
        name_var=strcat(name_var{1},name_var{2},name_var{3},name_var{4});
    else
        name_var = {name_f_save name_save_Var {'_binned'} {num2str(bin_size)} {'s.mat'}};
        name_var=strcat(name_var{1},name_var{2},name_var{3},name_var{4},name_var{5});

        num_exp = 2; % If the variable already exists, write _sampl_#
        while ~(exist(name_var{1}) == 0)

            name_var = {name_f_save name_save_Var {'_sampl_'} {num2str(num_exp)} {'_binned'} {num2str(bin_size)} {'s.mat'}};
            name_var=strcat(name_var{1},name_var{2},name_var{3},name_var{4},name_var{5},name_var{6},name_var{7});

            num_exp = num_exp + 1;

        end
    end

    % Save and free space
    save(name_var{1}, 'FSignal')
    clear fid
    clear name_var

end % function 