% FUNCTION HEADER begin
function [return_data] = exec_bm(setup_file_name, Process_run_text)
%if nargin < 3, repair_flag = false; end
disp('------------------started---------------------')
%maxNumCompThreads(8);
if ismac
    % Code to run on Mac plaform
elseif isunix
    % Code to run on Linux plaform
    GROUNDPATH = '~/QUENTIN/bm_program/';
elseif ispc
    % Code to run on Windows platform
 
    GROUNDPATH = 'C:/QUENTIN/bm_program/';
else
    disp('Platform not supported')
end

%setup_file = '/itp/MooseFS/dorn/run/site_1/setup_1_site_U1_t2_phsym.mat';
% if ~isdeployed()
%     addpath(GROUNDPATH)
%
%   addpath([GROUNDPATH, 'classes', SLASH])
%   addpath([GROUNDPATH, 'functions', SLASH])
%   addpath([GROUNDPATH, 'functions', SLASH, 'class_functions', SLASH])
%   addpath([GROUNDPATH, 'functions', SLASH, 'c_functions', SLASH])
%   addpath([GROUNDPATH, 'functions', SLASH, 'plot', SLASH])
% end

%% ============ LOAD PARAMETERS ==========================
%NOTE: include ./data/ in the setup_file_name


load(setup_file_name, 'parameters');
struct2vars(parameters)
%deploy Process_run number to ranges (two ranges), block number tells how many jobs a condor job should perform
cpt_active = false;


% ================ APPLY RUN ===============================================

Process_run = str2double(Process_run_text)+1;
if Process_run > numel(Vb_ran_cell)
    error(['Process run number to high! Choosen parameters need ', num2str(numel(Vg_ran_cell)), ' runs.']);
end
Vb_ran = Vb_ran_cell{Process_run};
Vg_ran = Vg_ran_cell{Process_run};
[Vbb, Vgg] = meshgrid(Vb_vec(:), Vg_vec(:));
% FUNCTION HEADER end


% =================SET UP=================


Geom = Geometry(basis.N_site, geometry.structure, geometry.bonds); % class
Bas = Basis(basis.N_site,basis.N_sector,basis.S_sector,basis.Docc_input,basis.N_spin); % class
if isfield(grid_mesh, 'uniform') && grid_mesh.uniform == false
    GF_grid = Grid(grid_mesh.points, grid_mesh.imag_i0plus);
else
    if ~isfield(grid_mesh, 'uniform')
        warning('add "uniform" option to grid_mesh in setup_file!')
    end
    GF_grid = Grid(linspace(grid_mesh.lim(1), grid_mesh.lim(2) ,grid_mesh.numel_points), grid_mesh.imag_i0plus); % class
end

if ~isfield(condor, 'method')
    warning('add "method" option to condor in setup_file!')
    condor.method = 'lindblad';
end


%Tab.s_print('E', {min(En.Energies)-1, energy_cut })

% ================BATH=================

if ~isfield(bath, 'numel_orbitals')
    warning('numel_orbitals field is missing in bath - old bath setup')
    bath.numel_orbitals = cellfun(@(x) size(x,1), bath.t_leads);
end
if ~isfield(bath, 'spin_conserved')
    warning('spin_conserved field is missing in bath - old bath setup')
    bath.spin_conserved = true(bath.N_leads,1);
end
if ~isfield(bath, 'spin_symmetric')
    warning('spin_symmetric field is missing in bath - old bath setup')
    bath.spin_symmetric = true(bath.N_leads,1);
end
if numel(bath.shift_bath) ~= bath.N_leads
    warning('shift bath in bath setup not defined for all baths - old bath setup')
    bath.shift_bath = repmat(bath.shift_bath,bath.N_leads,1);
end


bat = Bath(bath.GF, GF_grid, bath.shift_bath, bath.integrator, basis.N_spin, bath.numel_orbitals, bath.spin_conserved, bath.spin_symmetric);

%each run 1 output file!
II = zeros(1,1, basis.N_spin * bath.N_leads);
II_cpt = zeros(1,1, basis.N_spin * bath.N_leads);
II_bmcpt = zeros(1,1, basis.N_spin * bath.N_leads);
Trans_bmcpt = zeros(1,1, GF_grid.N);
Trans_cpt = zeros(1,1, GF_grid.N);

NN = zeros(1,1);
NN_full = cell(1,1);
Sigg = cell(1,1);
% NOTE: TO BE DELETED
%status_full = struct('no_solution',{},'multiple_solutions',{},'not_same_argument',{}, 'current_missmatch', {}, ...
%    'column_sum', {}, 'positive_eigenvalues', {}, 'eigenvalue_tolerance', {}, 'minimal_value', {}, 'energy_cut', {}, ...
%    'error_text', {}, 'choi_herm', {}, 'choi_neg_ew_ratio', {}, 'choi_mean_ew', {}) ;
%status_full(1,1).no_solution = [];


%TODO: delete, if all programs are changed - signalized with warning
if ~isfield(energy, 'eigenvalue_tolerance')
    energy.eigenvalue_tolerance = 10^-7;
    warning('energy.eigenvalue_tolerance not defined - adopt setup_file!!!')
end


%% ================== gate voltage ======================
for Vg_ind = 1:numel(Vg_ran)
    Vg = Vg_ran(Vg_ind);
    %TODO: a very special rule, develop generel rule for different geometries!
    if ~isfield(contact, 'gate_coupling') || isempty(contact.gate_coupling)
        system.b = system.b + eye(size(system.b)) * Vg;
    else
        system.b = system.b + diag(contact.gate_coupling*Vg);
    end
    % apply Vg to system.b
    
    if energy.fresh_calc || ~isfile(energy.file{Vg == Vg_vec})
        Hub = Hubbard(Geom, system.b, system.U, system.xi, system.V, Bas, system.part_hole_sym); % class
        tic
        En = Hub.f_solve(energy.calc, energy.second_cut);
        
        en_solve = toc;
        if condor.mes_flag
            disp(['Solving Energy: ', num2str(en_solve), 's, number of calculated eigenvalues: ', num2str(sum(En.Energies ~= 0)), ' of ', num2str(En.Numel_En)])
        end
        energy.energy_cut(Vg_vec == Vg) = min(En.Energies) + energy.first_cut;
        energy.energy_cut_2(Vg_vec == Vg) = min(En.Energies) + energy.second_cut;
        
        if Bas.N_spin == 2 && isempty(energy.spin_tag)
            %TODO: check out if energy is used!
            En = En.f_resort('NE');
        end
    else
        %NOTE: filename set up in setup_file_IDENTIFIER.m - may depend also on Energy cuts and blocksizes
        load(energy.file{Vg == Vg_vec},'En');  %dependent on Vg!!!!
    end
    
    if ~isempty(energy.block_energy_diff)
        En = condense_energy_blocks(En, energy.block_energy_diff);
        %NOTE: counts if the energy in an Average Energy Block is above the energy cut - if so, the energy cut is
        %adopted to match the gap boundary
        L = En.Energies >=En.Average_Energies - En.Gap_Energies & En.Energies<= En.Average_Energies+En.Gap_Energies & En.Average_Energies< energy.energy_cut(Vg_vec == Vg) & En.Energies >= energy.energy_cut(Vg_vec == Vg);
        
        if sum(L) > 0
            disp(['Energy cut raised since ', num2str(sum(L)), ' Energies were inside an Average Energy block but above the Energy cut'])
            energy.energy_cut(Vg_vec == Vg) = max(En.Energies(L))+ 10^-13;
        end
        Tab = Table(En, ['N',energy.spin_tag, 'A']); % class
        %sum(Tab.Block_Size.^2)
        %energy.energy_cut(Vg_vec == Vg) = En
        %energy.energy_cut_2(Vg_vec == Vg) =
        energy.energy_tag = 'A';
    else
        Tab = Table(En, ['N', energy.spin_tag, 'E']);
        energy.energy_tag = 'E';
    end
    
    
    % simple CPT
    %define an average temperature:
    %TODO: find a better way
    
    if cpt_active
        
        beta_system = mean(Beta(:));
        
        rho_diagonal = exp(-beta_system .*(En.Energies - min(En.Energies)));
        rho_diagonal(abs(rho_diagonal) < 10^-5) = 0;
        rho_diagonal = sparse(rho_diagonal./sum(rho_diagonal));
        
        block_index = Tab.Block_Index(energy.energy_tag, {-inf, energy.energy_cut(Vg_vec == Vg)});
        blocks = arrayfun(@(x) sparse(x,x), Tab.Block_Size, 'uniformoutput', false);
        block_size_cum = cumsum([Tab.Block_Size]);
        block_size_cum = [0; block_size_cum(1:end-1)];
        
        for k = 1:numel(block_index)
            block_id = block_index(k);
            diag_entries = rho_diagonal(block_size_cum(block_id) + (1:Tab.Block_Size(block_id)));
            if sum(diag_entries) > 10^-4
                blocks{block_id} = diag(diag_entries);
                %define unitary rotation
                if false
                    leng = Tab.Block_Size(block_id);
                    tr_to_have = trace(blocks{block_id});
                    
                    A = rand(leng,leng)-0.5 +  1i* rand(leng, leng);
                    A = A+A';
                    [U,V] = eig(A);
                    blocks{block_id} = U * blocks{block_id} * U';
                    blocks{block_id} = tr_to_have * 1/trace(blocks{block_id}) * blocks{block_id};
                end
            end
            
        end
        
        
        Sig_CPT = Sigma(En, blocks, Tab, energy.energy_cut(Vg_vec == Vg), energy.energy_cut_2(Vg_vec == Vg));
        spin_symmetric = true;
        [G_R_cpt, G_K_cpt] = Sig_CPT.f_greens_function_full(En, GF_grid, spin_symmetric, GF_grid.Imag_Infinitesimal);
        G_C_cpt = Green(GF_grid, G_R_cpt, G_K_cpt);
        %G_C_cpt.s_plot_spectral;
        %NOTE: takes ~90 sec
        %TODO: preinitialize with total Vb_vec or just relevant Vb_ran???? ->
        %TODO: add feature to reinitialize (extend the ranges of En and Vb
        
        %depends on Energy which changes!!!
    end
    
    % ----------------- Mu -----------------
    for Vb_ind = 1:numel(Vb_ran)
        
        Vb = Vb_ran(Vb_ind);
        disp([ 'Vb_ind: ', num2str(Vb_ind), ', Vb: ', num2str(Vb), ', Vg: ', num2str(Vg)])
        
        Energy_cut = energy.energy_cut(Vg_vec == Vg);
        Energy_cut_2 =  energy.energy_cut_2(Vg_vec == Vg);
        
        %Syntax for mu: Mu(Bath, Spin)
        %left bath minus Vb/2!
        %TODO: define a function how Vb is distributed to spins and baths
        Mu = [-Vb/2; Vb/2];
        bat = bat.f_preinitialize(En, Mu, Beta, energy.energy_cut(Vg_vec == Vg));
        bat.Mu = Mu;
        
        %En = En.f_resort('NES')
        %energy.energy_cut(Vg_vec == Vg) = energy.energy_cut(Vg_vec == Vg) +5;
        %energy.energy_cut_2(Vg_vec == Vg) = energy.energy_cut_2(Vg_vec == Vg) +5
        % condor.mes_flag = true
        run_id = [num2str(find(Vg == Vg_vec)), 'x', num2str(find(Vb == Vb_vec))];
        
        
        min_en_ind = energy.min_en_ind;
        
        K_file_name = ['./log_files/K_', condor.setup_id, '_', run_id, '.mat'];
        
        
        if exist(K_file_name, 'file') == 2
            load(K_file_name, 'K_data')
            struct2vars(K_data);
            % TODO: to be deleted
            K_recalc = true;
            
        end
        
        %[Sig, status] = get_sigma(K, Tab_cut, status, condor.mes_flag, Energy_cut, Energy_cut_2, 10^-10);
        %NOTE: really strange error!!!
        %[~, choi_herm, choi_ew] = Tab_cut.f_choi(K);
        %[K2, choi_herm] = Tab_cut.f_choi(K);
        %status.choi_herm = choi_herm;
        %status.choi_neg_ew_ratio = sum(real(choi_ew) < -10^-10)/sum(~isnan(choi_ew));
        %status.choi_mean_ew = mean(real(choi_ew(~isnan(choi_ew))));
        %disp(['number of nan elements: ', num2str(sum(isnan(choi_ew))/numel(choi_ew))])
        
        
        
        if min_en_ind == inf
            min_en_ind = Tab.Numel_Blocks;
        end
        if ~(exist(K_file_name, 'file') == 2) || ((exist('K_recalc','var') == 1) && K_recalc == true)
            
            switch condor.method
                case 'bm'
                    symmetrized_version = false;
                    [K, I, Sig, bat, Tab_cut, status, calc_time, column_sum, ev_min, choi, X_cur] = born_markov_full_dynamic(En, contact.contact_site, Tab, Mu, contact.coupling, Beta, bat, condor, Energy_cut, Energy_cut_2, run_id, min_en_ind, energy.eigenvalue_tolerance, symmetrized_version);
                case 'lindblad'
                    multi_symmetrization = false;
                    [K, I, Sig, bat, Tab_cut, status, calc_time, column_sum, ev_min, choi, X_cur, H_LS_calculated] = lindblad(En, contact.contact_site, Tab, Mu, contact.coupling, Beta, bat, condor, Energy_cut, Energy_cut_2, run_id, min_en_ind, energy.eigenvalue_tolerance, multi_symmetrization);
                case 'lindblad_multi'
                    multi_symmetrization = true;
                    [K, I, Sig, bat, Tab_cut, status, calc_time, column_sum, ev_min, choi, X_cur, H_LS_calculated] = lindblad(En, contact.contact_site, Tab, Mu, contact.coupling, Beta, bat, condor, Energy_cut, Energy_cut_2, run_id, min_en_ind, energy.eigenvalue_tolerance, multi_symmetrization);
                    
                    
                case 'bm_sym'
                    symmetrized_version = true;
                    [K, I, Sig, bat, Tab_cut, status, calc_time, column_sum, ev_min, choi, X_cur] = born_markov_full_dynamic(En, contact.contact_site, Tab, Mu, contact.coupling, Beta, bat, condor, Energy_cut, Energy_cut_2, run_id, min_en_ind, energy.eigenvalue_tolerance, symmetrized_version);
                    
            end
            %{
            symmetrized_version = false;
                    [K, I, Sig, bat, Tab_cut, status, calc_time, column_sum, ev_min, choi, X_cur] = born_markov_full_dynamic(En, contact.contact_site, Tab, Mu, contact.coupling, Beta, bat, condor, Energy_cut, Energy_cut_2, run_id, min_en_ind, energy.eigenvalue_tolerance, symmetrized_version);
             multi_symmetrization = false;
                    [K2, I2, Sig2, bat, Tab_cut, status2, calc_time2, column_sum, ev_min, choi, X_cur] = lindblad(En, contact.contact_site, Tab, Mu, contact.coupling, Beta, bat, condor, Energy_cut, Energy_cut_2, run_id, min_en_ind, energy.eigenvalue_tolerance, multi_symmetrization);
              symmetrized_version = true;
                    [K3, I3, Sig3, bat, Tab_cut, status3, calc_time3, column_sum, ev_min, choi, X_cur] = born_markov_full_dynamic(En, contact.contact_site, Tab, Mu, contact.coupling, Beta, bat, condor, Energy_cut, Energy_cut_2, run_id, min_en_ind, energy.eigenvalue_tolerance, symmetrized_version);
                      
            
            
            
            %}
            
            
            
            %[K, I, Sig, bat, Tab_cut, status, calc_time, column_sum, ev_min, choi, Sig_result] = lindblad(En, contact.contact_site, Tab, Mu, contact.coupling, Beta, bat, condor, Energy_cut, Energy_cut_2, run_id, min_en_ind, energy.eigenvalue_tolerance);
            
            %[K2, I2, Sig2, bat, Tab_cut, status] = born_markov_full_dynamic(En, contact.contact_site, Tab, Mu, contact.coupling, Beta, bat, condor, Energy_cut, Energy_cut_2, run_id, min_en_ind, energy.eigenvalue_tolerance);
            
            %             disp(['Multiple solutions: ', num2str(status.multiple_solutions), ', minimal eigenvalue: ', num2str(status.minimal_value)])
            %             figure
            %             ew = eig(full(K)); plot(real(ew), imag(ew), 'xk'), grid on
            %             status
            Sig.s_dom_blocks
            %analyse_K(K, Tab_cut, calc_time, Sig, status, Energy_cut, Energy_cut_2)
            
            Tab_cut.Energy = [];
            if isfield(condor, 'save_final_K') && condor.save_final_K
                
                K_recalc = false;
                
                
                K_data = vars2struct('K', 'I', 'Tab_cut',  'Sig', 'status', 'Energy_cut', 'Energy_cut_2', 'calc_time', 'Vg', 'Vb', 'Vb_vec', 'Vg_vec', 'setup_file_name', 'K_recalc', 'X_cur');
                if isfield(condor, 'return_data') && condor.return_data
                    return_data.KK{Vg_ind, Vb_ind} = K;
                    return_data.Sig{Vg_ind,Vb_ind} = Sig;
                    return_data.Tab_cut{Vg_ind, Vb_ind} = Tab_cut;
                    return_data.Energy_cut(Vg_ind, Vb_ind) = Energy_cut;
                    return_data.Energy_cut_2(Vg_ind, Vb_ind) = Energy_cut_2;
                    return_data.calc_time{Vg_ind, Vb_ind} = calc_time;
                    return_data.Vg(Vg_ind, Vb_ind) = Vg;
                    return_data.Vb(Vg_ind, Vb_ind) = Vb;
                    return_data.Vb_vec = Vb_vec;
                    return_data.Vg_vec = Vg_vec;
                    
                else
                    save(K_file_name, 'K_data')
                end
            end
            
        end
        
        if status.no_solution
            disp('===================NO SOLUTION!!!!===================')
        end
        disp('================ born_markov_full_dynamic finished ===============')
        %[K, I, Sig, bat, status] = born_markov_full(En, contact.contact_site, Tab, Mu, contact.coupling, Beta, bat, condor.mes_flag, energy.energy_cut(Vg_vec == Vg), energy.energy_cut_2(Vg_vec == Vg));
        
        %add condor.save_final_K, condor.save_intermediate_K
        
        %end
        
        % diag_index = Tab.Diagonal_Index();
        % K_relevant = K(diag_index,:);
        
        %[K2,I,Sig] =  born_markov_2(En, contact.contact_site, Tab, Mu, repmat([contact.coupling{1}(1); contact.coupling{2}(2)],1,2), Beta(1), bat, condor.mes_flag, energy.energy_cut(Vg_vec == Vg));
        
        weight = Sig.s_dom_blocks;
        Sig.Tab.Energy = []; %delete Energy from Sig, to save storage.
        
        
        NN = sum(weight.N_part.*weight.ev_weight);
        NN_full =  Sig.f_diagonal;
        Sigg = Sig;
        II = I(:);
        status_full = status;
        
        
        % CPT
        
        %Sig.f_greens_function_full(grid, spin_symmetric, i0p);
        if ~ismember('S', En.Structure)
            spin_symmetric = false;
        else
            % TODO: implement spin_symmetric as configuration parameter! derive from Bath class
            spin_symmetric = true;
        end
        %all calculated H_LS blocks
        %H_LS_calculated( Sig.Tab.Old_Indices)
        
        if cpt_active
            
            [G_R, G_K] = Sig.f_greens_function_full(En, GF_grid, spin_symmetric, GF_grid.Imag_Infinitesimal); %TODO: make it bigger than mesh?
            
            %G_C = Green(grid, G_R(1:8,1:8,:), G_K(1:8,1:8,:));
            G_C = Green(GF_grid, G_R, G_K);
            
            [Current_bmcpt, Detail_Current, I_bmcpt, FL_detail, Transmission_bmcpt] = CPT_current(G_C, bat, contact.coupling, GF_grid, Mu, Beta, Bas, En);
            [Current_cpt, Detail_Current, I_cpt, FL_detail, Transmission_cpt] = CPT_current(G_C_cpt, bat, contact.coupling, GF_grid, Mu, Beta, Bas, En);
            
            II_bmcpt  = reshape(I_bmcpt, bat.Numel_Baths*Bas.N_spin,1);
            II_cpt = reshape(I_cpt, bat.Numel_Baths * Bas.N_spin,1);
            Trans_bmcpt = Transmission_bmcpt(:);
            Trans_cpt = Transmission_cpt(:);
            
        else
            I_bmcpt = nan(size(I));
            I_cpt = nan(size(I));
            
        end
        
        disp('==============================Current: ====================================')
        disp([ 'Vb_ind: ', num2str(Vb_ind), ', Vb: ', num2str(Vb), ', Vg: ', num2str(Vg)])
        table( -I(1:2,1), -I_bmcpt(1:2,1), -I_cpt(1:2,1), 'variablenames', {'BM', 'BMCPT', 'CPT'})
        
        
        %NOTE: condor.file_name is composed of data_IDENTIFIER_date_setup_version
        %NOTE: switch to internal run number instead of condor run number (Process_run_text)
        
        internal_run_number = find(Vb == Vbb(:) & Vg == Vgg(:));
        file_name = [condor.file_name, condor.file_name_pre_number, num2str(internal_run_number) , '.mat'];
        %file_name = [condor.file_name, condor.file_name_pre_number, Process_run_text , '.mat'];
        
        data = vars2struct('parameters', 'Vb_ran', 'Vg_ran', 'Vb', 'Vg', 'II', 'NN', 'NN_full', 'II_cpt', 'II_bmcpt', ...
            'status_full', 'Sigg', 'Trans_cpt', 'Trans_bmcpt');
        if (isfield(condor,'return_data') && condor.return_data  )
            return_data.data = data;
        else
            save(file_name, 'data')
            return_data = [];
            disp(['finished, file saved: ' file_name])
        end
        
    end
end


% figure,
%
% plot(Vb_ran, II(:,:,2))
% hold on
% plot(Vb_ran, II_cpt(:,:,2))
% plot(Vb_ran, II_bmcpt(:,:,2))



disp(['Condor / Slurm run finished: ', Process_run_text])



end
