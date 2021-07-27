
function post_run_fuse_files(ID, save_file_name)


% how to deal with not finished condor / slurm jobs:
% a) always create a new combined file (even though the newer contains the old data)
% b) if some jobs were not completed, reshuffle condor and print the number of runs to be done 
% c) if all jobs were completed, reset the condor jobs to the full amount

% Future aim: start post_run_fuse_files automatically after all Condor jobs are done and print the result in a
% special folder (one status file for each setup_file, for each Condor run a line)
% TODO: think about situation when more jobs are started than needed, 



if nargin < 2 || isempty(save_file_name)
    date_format = 'yyyy-mm-dd';
    time_format = 'HH:MM';
    date_string = datestr(now, [date_format, '_', time_format]);
    
    save_file_name = ['./data/combined_data_', ID, '_', date_string, '.mat'];
end






% load setup file
setup_file_name = ['./data/setup_', ID, '.mat'];
load(setup_file_name, 'parameters')
struct2vars(parameters)


%NOTE: find out which internal run ids should be there, or check by listings!
total_number_of_runs = numel(Vb_vec) * numel(Vg_vec);
[VVB, VVG] = meshgrid(Vb_vec, Vg_vec);
internal_run_id = cell(numel(Vb_ran_cell),1);
for u = 1:numel(Vb_ran_cell)
    [test_Vb, test_Vg] = meshgrid(Vb_ran_cell{u}, Vg_ran_cell{u});
    internal_run_id{u} = arrayfun(@(x,y) find(x == VVB & y == VVG), test_Vb, test_Vg);
    
end

int_run_id = [internal_run_id{:}]; %as list


% NOTE: find out which runs have been finished (file ids)
listing = dir(['./data/data_', ID, '_*.mat']);
if isempty(listing)
    actual_file_id = [];
else
    [actual_file_id, ~] = regexp([listing(:).name], ['data_', ID, '_(\d*).mat'], 'tokens', 'match');
    actual_file_id = ([actual_file_id{:}]);
    actual_file_id = cellfun(@(x) str2double(x), actual_file_id(:));

end
%TODO: old to be deleted, just for convergence reasons - check with setup time 190906
if ismember(0, actual_file_id)
    warning('old data')
    actual_file_id = actual_file_id +1;
end

missing_file_id = 1:total_number_of_runs;
missing_file_id(actual_file_id) = [];


if grid_mesh.uniform
    grid_N = grid_mesh.numel_points;
else
    grid_N = numel(grid_mesh.points);
end

% NOTE: Initialize
% check if there is already calculated data (OLD - to be deleted in future releases)

file_name_preprocessed = [condor.file_name(1:7), 'combined_', condor.file_name(8:end), '-1.mat'];
if exist(file_name_preprocessed, 'file')
    load(file_name_preprocessed, 'data_combined')
    struct2vars(data_combined)
    % struct2vars(parameters)
    %TODO: check if dimensions fit - if parameters coincide
else
    II = nan(numel(Vg_vec), numel(Vb_vec), basis.N_spin * bath.N_leads);
    II_cpt = nan(numel(Vg_vec), numel(Vb_vec), basis.N_spin * bath.N_leads);
    II_bmcpt = nan(numel(Vg_vec), numel(Vb_vec), basis.N_spin * bath.N_leads);
    VB = nan(numel(Vg_vec), numel(Vb_vec));
    VG = nan(numel(Vg_vec), numel(Vb_vec));
    NN = nan(numel(Vg_vec), numel(Vb_vec));
    NN_full = cell(numel(Vg_vec), numel(Vb_vec));
    Sigg = cell(numel(Vg_vec), numel(Vb_vec));
    Trans_cpt = nan(numel(Vg_vec), numel(Vb_vec), grid_N);
    Trans_bmcpt = nan(numel(Vg_vec), numel(Vb_vec), grid_N);
    %TODO: change to status_full
    no_solution = nan(numel(Vg_vec), numel(Vb_vec));
    multiple_solutions = nan(numel(Vg_vec), numel(Vb_vec));
    not_same_argument =  nan(numel(Vg_vec), numel(Vb_vec));
    current_missmatch = nan(numel(Vg_vec), numel(Vb_vec));
    column_sum = nan(numel(Vg_vec), numel(Vb_vec));
    positive_eigenvalues = nan(numel(Vg_vec), numel(Vb_vec));
    eigenvalue_tolerance = nan(numel(Vg_vec), numel(Vb_vec));
    minimal_value = nan(numel(Vg_vec), numel(Vb_vec));
    energy_cut = nan(numel(Vg_vec), numel(Vb_vec));
    partial_trace_converged = nan(numel(Vg_vec), numel(Vb_vec));
    energy_domain_limit = nan(numel(Vg_vec), numel(Vb_vec));
    steady_state_negative = nan(numel(Vg_vec), numel(Vb_vec));
    choi_herm = nan(numel(Vg_vec), numel(Vb_vec));
    nu_m = nan(numel(Vg_vec), numel(Vb_vec));
    nu_t = nan(numel(Vg_vec), numel(Vb_vec));
    full_gamma_positivity = nan(numel(Vg_vec), numel(Vb_vec))
    
    
    error_text = cell(numel(Vg_vec), numel(Vb_vec));
    
    
    
end

for k = 1:numel(actual_file_id)
    disp(['file #', num2str(k), ', file ID: ', num2str(actual_file_id(k))])

    file_name = [condor.file_name, condor.file_name_pre_number, num2str(actual_file_id(k)), '.mat'];
    if exist(file_name, 'file') == 2
        load(file_name, 'data')
        
        if ~isfield(data, 'Vb')
            data.Vb = data.Vb_ran;
        end
        if ~isfield(data, 'Vg')
            data.Vg = data.Vg_ran;
        end
       
               
        [x_control,ix] = ismember(data.Vg, Vg_vec);
        if any(~x_control)
            % merge!!! insert new values to Vg_vec, how could this happen? maybe do extra with extra merge
            error('to be done')
        end
        [y_control,iy] = ismember(data.Vb,Vb_vec);
        if any(~y_control)
            % merge!!! insert new values to Vb_vec, how could this happen? maybe do extra with extra merge
            error('to be done')
        end
        
        II(ix,iy,:) = data.II;
        II_cpt(ix,iy,:) = data.II_cpt;
        II_bmcpt(ix,iy,:) = data.II_bmcpt;
        NN(ix,iy,:) = data.NN;
        NN_full{ix,iy} = data.NN_full;
        if isfield(data, 'Trans_cpt')
            Trans_cpt(ix,iy,:) = data.Trans_cpt;
        end
        if isfield(data,'Trans_bmcpt')
            Trans_bmcpt(ix,iy,:) = data.Trans_bmcpt;
        end

        no_solution(ix,iy) = data.status_full.no_solution;
        multiple_solutions(ix,iy) = data.status_full.multiple_solutions;
        not_same_argument(ix,iy) = data.status_full.not_same_argument;
        current_missmatch(ix,iy) = data.status_full.current_missmatch;
        column_sum(ix,iy) = data.status_full.column_sum;
        positive_eigenvalues(ix,iy) = data.status_full.positive_eigenvalues;
        eigenvalue_tolerance(ix,iy) = data.status_full.eigenvalue_tolerance;
        minimal_value(ix,iy) = data.status_full.minimal_value;
        energy_cut(ix,iy) = data.status_full.energy_cut;
        if ~isempty(data.status_full.partial_trace_converged)
            partial_trace_converged(ix,iy) = data.status_full.partial_trace_converged;
        end
        if ~isempty(data.status_full.energy_domain_limit)
            energy_domain_limit(ix,iy) =  data.status_full.energy_domain_limit;
        end
        steady_state_negative(ix,iy) = data.status_full.steady_state_negative; 
        choi_herm(ix,iy) = data.status_full.choi_herm;
        nu_m(ix,iy) = data.status_full.nu_m;
        nu_t(ix,iy) = data.status_full.nu_t;
        full_gamma_positivity(ix,iy) = data.status_full.full_gamma_positivity;
    
        for l = 1:numel(iy)
            error_text{ix,iy(l)} = data.status_full(l).error_text;
        end
        Sigg{ix,iy} = data.Sigg;
        VB(ix,iy) = data.Vb;
        VG(ix,iy) = data.Vg;
        


        
        %console_file = ['./log_files/Console', num2str(k), '.out'];
        %fileID = fopen(console_file, 'r');
        %text = fscanf(fileID, '%c',[1,inf]);
        %fclose(fileID);
        %warn.multiple_solutions(ix,iy) = numel(regexp(text,'WARNING: multiple solutions, take smallest absolute value'));
        %warn.no_solution(ix,iy) = numel(regexp(text, 'WARNING: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NO SOLUTION!!!!!!!!!!!!!!!!!!!!!!!!')); 
        
    else
        warning('existing file could not be opened - check it!')
    end
end




data_combined = vars2struct('II','II_cpt', 'II_bmcpt', 'NN', 'NN_full', 'VB', 'VG', 'Sigg', 'parameters', ...
    'no_solution', 'multiple_solutions', 'not_same_argument', 'current_missmatch', 'column_sum', 'positive_eigenvalues', ...
    'minimal_value', 'eigenvalue_tolerance', 'energy_cut', 'partial_trace_converged', 'energy_domain_limit', 'error_text', ...
    'steady_state_negative',  'choi_herm', 'nu_m', 'nu_t', 'full_gamma_positivity', 'Trans_cpt', 'Trans_bmcpt');


%check which runs are missing:
if sum(isnan(VB(:))) > 0 %indicates that overall some runs as planned in setup were not successful
    disp([num2str(numel(find(isnan(VB)))), ' not finished: ', num2str(find(isnan(VB(:))).')])

    % reshuffle condor
    L = cellfun(@(x) any(ismember(x, missing_file_id)), internal_run_id);
    Vb_ran_cell(~L) = [];
    Vg_ran_cell(~L) = [];
    disp(['Number of runs to perform: ', num2str(numel(Vb_ran_cell))]);
    parameters = vars2struct('condor', 'Vg_vec', 'Vb_vec', 'geometry', 'basis',...
                         'grid_mesh', 'energy', 'system', 'contact', 'bath', 'Beta', 'Vg_ran_cell', 'Vb_ran_cell');
  
    save([setup_file_name],'parameters')
    
else
    if ~all(ismember(1:total_number_of_runs, int_run_id))
        %reset Vb_ran_cell and Vg_ran_cell
        if contact.couple_bias_gate
            Vg_vec = Vb_vec;
            Vg_block_size = condor.block_size;
            number_full_Vg_blocks = floor(numel(Vg_vec)/Vg_block_size);
            Vg_block_size_vector = [ones(1,number_full_Vg_blocks)*Vg_block_size, numel(Vb_vec) - number_full_Vg_blocks * Vg_block_size];
            Vb_ran_cell = mat2cell(Vb_vec,1, Vg_block_size_vector );
            Vg_ran_cell = mat2cell(Vg_vec,1, Vg_block_size_vector );
        else
            if condor.block_size > numel(Vb_vec)
                condor.block_size = numel(Vb_vec) * floor(condor.block_size / numel(Vb_vec));
                Vg_block_size = condor.block_size/numel(Vb_vec);
            else
                Vg_block_size = 1;
            end
            number_full_Vg_blocks = floor(numel(Vg_vec)/Vg_block_size); %without the rest
            Vg_block_size_vector = [ones(number_full_Vg_blocks,1)*Vg_block_size; ones(mod(numel(Vg_vec),Vg_block_size)>0,1)];
            if ~isdeployed(), disp(['Numel of Vg_blocks in total: ', num2str(numel(Vg_block_size_vector))]); end
            
            number_Vb_blocks = floor(numel(Vb_vec) /condor.block_size); % without the rest
            last_block_size_Vb = numel(Vb_vec) - condor.block_size*number_Vb_blocks;
            Vb_block_size_vector = [ones(1,number_Vb_blocks)*condor.block_size,last_block_size_Vb*ones(1, last_block_size_Vb > 0)];
            if ~isdeployed(), disp(['Numel of Vb_blocks per Vg_entry: ', num2str(numel(Vb_block_size_vector))]); end
            
            Vb_ran_cell = mat2cell(repmat(Vb_vec, numel(Vg_block_size_vector),1),ones(size(Vg_block_size_vector)),Vb_block_size_vector);
            Vg_ran_cell = mat2cell(repmat(Vg_vec(:), 1, numel(Vb_block_size_vector)), Vg_block_size_vector, ones(size(Vb_block_size_vector)));
        end
        
        parameters = vars2struct('condor', 'Vg_vec', 'Vb_vec', 'geometry', 'basis',...
                         'grid_mesh', 'energy', 'system', 'contact', 'bath', 'Beta', 'Vg_ran_cell', 'Vb_ran_cell');
  
        save(setup_file_name,'parameters')
    end
end
    
%check which condor / slurm runs failed

    

%listing = dir(['./data/combined_data_benzene', '_', date_format, '*.mat']);



save(save_file_name, 'data_combined')
disp(['file saved: ', save_file_name])
%load(save_file_name, 'data_combined')



end

