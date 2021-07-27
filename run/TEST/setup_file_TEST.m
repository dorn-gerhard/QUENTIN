%SET UP FILE
% =================SET UP=================

function [setup_id] =  setup_file_TEST(setup)


IDENTIFIER = 'TEST';

% =============================== Basis ===================================
basis.N_spin = 2;
basis.N_site = 3;
basis.N_sector = []; %empty for all
basis.S_sector = []; %empty for all
basis.Docc_input = 1; % 1 for all, 0 for up and down, 2 for none, up and down; 3 for up, down and double
 % Table for code of Docc_input
            %         |    0    |    1    |    2    |
            % ---------------------------------------
            %     1   |    y    |    y    |    y    |
            %     0   |    n    |    y    |    n    |
            %     2   |    y    |    y    |    n    |
            %     3   |    n    |    y    |    y    |
            %         |    y    |    n    |    y    |

% ============================== System ===================================

system.U = [3 0 0; 0 3 0; 0 0 3];
system.V = [0 0 0 ; 0 0 0 ; 0 0 0];
system.b = [0 -1 -1; -1 0.01 -1; -1 -1 0];
system.xi = []; %not using  xi

system.part_hole_sym = true;


%Hub_spectrum(basis,system)
% ============================== Energy ===================================

%energy calc used to calculate the eigenvalues of the system
energy.calc.tolerance = 10^-14; %used for rounding eigenvalues and determining spin and particle sectors
energy.calc.max_number = 50; % a set energy.second_cut overrules the max_number!

energy.eigenvalue_tolerance = 10^-7;

energy.first_cut = inf;
energy.second_cut = inf;

energy.fresh_calc = false; % since we use files

energy.min_en_ind = 100; 
energy.block_energy_diff =  10^-10; %empty for not adding up block energies
energy.spin_tag = 'S';  
%NOTE: blocks to take into account, if inf, full K is established, if empty standard value is taken, determines the minimal size of K



% ============================= Geometry ==================================
geometry.structure = 'f'; %'l' ... linear, 'c' ... cyclic, 'f' ... full
geometry.bonds = ones(basis.N_site);

% ============================= Condor ====================================
condor.block_size = 1; % for condor
condor.mes_flag = true;

condor.save_final_K = true;
condor.save_intermediate_K = false;
condor.method = 'lindblad'; %'bm'; %'lindblad_multi'; % 'bm_sym'; %
if isfield(setup, 'research_question' )
    condor.research_question = "How is the current characteristics?";
else
    condor.research_question = '';
end
condor.research_answer = '';

condor.increase_tolerance = 0; % number of blocks from break set (steps of 10) after which tolerance is increased by one order (*10)

% ============================ Gate Voltage ===============================
Vg_delta = 0.5;
Vg_vec = 0; %0:Vg_delta:10;

% ============================ Bias Voltage ===============================
Vb_delta = 0.1;
Vb_vec = 0:Vb_delta:3;


% ============================== Grid =====================================
grid_mesh.uniform = false;
grid_mesh.lim = [-6,6];
grid_mesh.numel_points = 24001; %24001 
grid_mesh.imag_i0plus = 10^-7;
grid_mesh.points = [-20:-10, -9.5:0.5:-6.5, linspace(-6,6,19345), 6.5:0.5:9.5, 10:20];


% ============================= Baths =====================================
%To be defined: bath.GF, bath.shift_bath, bath.integrator, bath.numel_orbitals, bath.spin_conserved, bath.spin_symmetric)

% either assume wide band limit, thight binding or a Green's function from some input file
if grid_mesh.uniform
    bath.type = 'from_file';
    GF_input = {G_BL_ret; G_BR_ret}; %''; %from file
    bath.GF = GF_input;
    bath.numel_orbitals = cellfun(@(x) size(x,1), bath.GF);
    
else
    bath.type = 'tight_binding';%'wide_band'; %'from_file'; % 
    bath.GF = {};
    % ----------------------two orbitals (numel_orbitals = 2), different in both baths -------------
    % ----------------------for spin dependent baths - set lead_coupling to 0 and numel_orbitals to 1
%     lead_coupling = 0.1;    %left                        %right
%     bath.t_leads =            {[0.25, lead_coupling   ;...
%                           conj(lead_coupling), 1.25]; ...
%                                                     [0.25, lead_coupling; ...
%                                                      conj(lead_coupling), 1.25]};
%     bath.eps_0 =              {[0,   0; ...
%                            0  ,  -2.15]; ...
%                                                     [0, 0; ...
%                                                      0, 1.15]};
    
    
    % ------------------------ only one orbital - different for both baths
%     t_leads =  {6; 4};
%     eps_0 = {0;2};
%     
%     bath.t_leads = t_leads;
%     bath.eps_0 = eps_0;
%     
    
    % ------------------------ only one orbital - the same both baths
    t_leads =  {6; 6};
    eps_0 = {0;0};
    
    bath.t_leads = t_leads;
    bath.eps_0 = eps_0;
    
    
    % ------------------------ define numel orbitals
    bath.numel_orbitals = cellfun(@(x) size(x,1), bath.t_leads);
    
    
end

bath.integrator = 'Fast_Detail'; %'Fast', 'Matlab', 'Detail', 'Constant'

% for wide band and tight binding leads
bath.N_leads = 2;

bath.spin_conserved = [true; true]; %mixture of spins (mostly true) - make an example for spin_conserved = false

bath.spin_symmetric = [true; true]; %same bath for spin up and down (spin_conserved necessary)

bath.shift_bath =     [false; false];




% bath.t_leads = cell(bath.N_leads, 1); %rows: baths, columns: spin
% bath.eps_0   = cell(bath.N_leads,  1); %rows: baths, columns: spin
% [bath.t_leads{:}] = deal(t_leads); %same bath on both sides
% [bath.eps_0{:}] = deal(eps_0);



bath.wideband.cutoff = true; % if true, real G_ret is created
bath.wideband.reconstruct_real_part = true; % if true real part is calculated by Kramers Kronig

%   bath_temp = Bath_New.c_wide_band(GF_grid, bath.t_leads, bath.shift_bath, bath.eps_0, bath.wideband.cutoff, bath.wideband.reconstruct_real_part, bath.integrator);
%   bath_temp = Bath_New.c_tight_binding(GF_grid, bath.t_leads, bath.shift_bath, bath.eps_0, bath.integrator);




% ============================= Contact ===================================
%Syntax for coupling: Coupling{Bath, Spin}(System Site, Bath orbital) !!!!

contact.couple_bias_gate = false;


coupling_factor = 0.1;

%NOTE: if spin is equal in bath, don't use it (as for Bath Green function)
if bath.numel_orbitals(1) == 2
    contact.coupling{1,1} = [1, 0.4;  0, 0;  0, 0] * coupling_factor;
else
    contact.coupling{1,1} = [1; 0; 0] * coupling_factor;
end


if bath.numel_orbitals(2) == 2
    contact.coupling{2,1} = [0, 0;  0, 0;  0.5, 1] * coupling_factor;
else
    contact.coupling{2,1} = [0; 0; 1] * coupling_factor;
end

        
        
contact.contact_site = 1:basis.N_site;
contact.gate_coupling = ones(1,basis.N_site);

% ============================= Temperature ===============================

beta = 50;
Beta = [beta,beta;beta,beta];% Syntax for Beta: Beta(Bath, [Spin - not meaningful])



% ============================= Done with setup ===========================

if ~isdeployed()
	%GROUNDPATH = '~/bm_program/';
	%GROUNDPATH = './bm_program';
    GROUNDPATH = '/QUENTIN/bm_program/';

	addpath(GROUNDPATH)
	addpath([GROUNDPATH, 'classes/'])
	addpath([GROUNDPATH, 'functions/'])
	addpath([GROUNDPATH, 'functions/class_functions/'])
	addpath([GROUNDPATH, 'functions/c_functions/'])
	addpath([GROUNDPATH, 'functions/plot/'])
end


% ======================== create directories =============================
% create directories ./data/ and ./log_files/ if not already there

if ~(exist('data', 'dir') == 7)
    mkdir('data')
end


if ~(exist('log_files', 'dir') == 7)
    mkdir('log_files')
end

% ===================== Setup Grid and Bath ===============================

if  grid_mesh.uniform == false
    GF_grid = Grid(grid_mesh.points, grid_mesh.imag_i0plus);
else
    GF_grid = Grid(linspace(grid_mesh.lim(1), grid_mesh.lim(2) ,grid_mesh.numel_points), grid_mesh.imag_i0plus); % class
end

Bath.f_check_Coupling(bath, contact.coupling, basis)



if strcmp(bath.type,'from_file')
    
elseif strcmp(bath.type, 'tight_binding')
    bath_temp = Bath.c_tight_binding(GF_grid, bath.t_leads, bath.shift_bath, bath.eps_0, bath.integrator, basis.N_spin, bath.numel_orbitals);
    bath.GF = bath_temp.GF;
elseif strcmp(bath.type, 'wide_band')
    bath_temp = Bath.c_wide_band(GF_grid, bath.t_leads, bath.shift_bath, bath.eps_0, bath.wideband.cutoff, bath.wideband.reconstruct_real_part, bath.integrator);
    bath.GF = bath_temp.GF;
end




% plot_GF(bath.GF{1}, GF_grid)
%plot_GF(bath_temp.GF{1}, GF_grid)
            
%g = Green(GF_grid, bath_temp.GF{1});
% create GF from wide_band



% if lower(bath.type(1)) == 't' %tight binding
%     bat = Bath_New.c_tight_binding(grid, bath.t_leads, bath.shift_bath);
% elseif lower(bath.type(1)) == 'w' %wideband limit
%     bat = Bath_New.c_wide_band(grid, bath.t_leads, bath.shift_bath, ...
%         bath.eps_0, bath.wideband.cutoff, bath.wideband.reconstruct_real_part);
% elseif lower(bath.type(1)) == 'f' %from file
%     bat = Bath_New(bath.GF, grid, bath.shift_bath);
% end
%     


% -------------------------- Setup Condor ---------------------------------

date_format = 'yymmdd';
%time_format = 'HH:MM';  % optional
date_string = datestr(now, [date_format]);

%check if setup_file already exists:
listing = dir(['./data/setup_', IDENTIFIER, '_', date_string, '_*.mat']);
number = numel(listing);

setup_id = [IDENTIFIER, '_', date_string, '_', num2str(number)];

save_file_name = ['./data/combined_data_',IDENTIFIER, '_', date_string, '_', num2str(number), '.mat'];



condor.file_name = ['./data/data_', setup_id];
condor.file_name_pre_number = '_';
condor.setup_id = setup_id;

% setup condor ranges
%NOTE: addressing scheme is always in such a way, so you can fill in a square manner:
%II(x, y, :) =  ... so all ways corresponding to any combination of x and y are calculated.
% (1,1), (1,2), (1,3), (2,1) is not possible -> instead just the first row is chosen.
% =================================================================
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

if ~isdeployed(), disp(['Numel of runs: ', num2str(numel(Vb_ran_cell))]); end
% =====================================================================


%load 'PtBDTPt_channelI_LDA+U.mat'            
%Geom = Geometry(basis.N_site, geometry.structure); % class
%Bas = Basis(basis.N_site,basis.N_sector,basis.S_sector,basis.Docc_input,basis.N_spin); % class
%addpath '/afs/itp.tugraz.at/user/dorn/Matlab/CPT/classes/Benzene' %not necessary any more

%load 'Hamiltonian_data.mat'
%load 'PtBDTPt_channelI_LDA.mat'
%load 'PtBDTPt_channelI_Ufull_V2.0.mat'
%load 'PtBDTPt_channelI_Ufull.mat' % OLD VERSION



% ------------------------------- Setup Energy ----------------------------

%grid = Grid(linspace(grid_mesh.lim(1), grid_mesh.lim(2) ,grid_mesh.numel_points), grid_mesh.imag_i0plus); % class


%TODO: think about calculating the energy_cuts dynamically if energy is fresh calculated!
% set just in setup file
energy.energy_cut = inf(1, numel(Vg_vec));
energy.energy_cut_2 = inf(1, numel(Vg_vec));

if ~energy.fresh_calc
    precalculate_now = true;
    save_energy_now = true;
    
    
    Bas = Basis(basis.N_site,basis.N_sector,basis.S_sector,basis.Docc_input,basis.N_spin); % class
    Geom = Geometry(basis.N_site, geometry.structure, ones(basis.N_site)); % class
    
    
    n_groundstates = min(1000, Bas.Numel_Bas);
    Groundstate = zeros(numel(Vg_vec), n_groundstates);
    for k = 1:1:numel(Vg_vec)
        Vg = Vg_vec(k);
        energy.file{k} = ['./data/energy_', setup_id '_Vg_', num2str(Vg), '.mat'];
        
        
        if precalculate_now
            Hub = Hubbard(Geom, system.b + contact.gate_coupling*Vg, system.U, system.xi, system.V, Bas, system.part_hole_sym); % class
            tic
            
            En = Hub.f_solve(energy.calc, energy.second_cut); % class
            %Q-matrix preinitialize
            
            en_solve = toc;
            disp(['Solving Energy: ', num2str(en_solve), 's, number of calculated eigenvalues: ', num2str(sum(En.Energies ~= 0)), ' of ', num2str(En.Numel_En)])

            
            EE = sort(En.Energies);
            Groundstate(k,:) = EE(1:n_groundstates);
            disp(['Vg: ', num2str(Vg_vec(k)), ', Groundstate: ', num2str(min(En.Energies)), ', first energy cut: ', ...
                num2str(sum(EE < (EE(1) + energy.first_cut))), ' states, second energy cut: ', num2str(sum(EE < (EE(1) + energy.second_cut))), ' states'])
            energy.energy_cut(k) = Groundstate(k,1) + energy.first_cut;
            energy.energy_cut_2(k) = Groundstate(k,1) + energy.second_cut;
            
            En = En.f_Qmatrix_generate( energy.energy_cut_2(k));
            
            if save_energy_now
                save(energy.file{k}, 'En')
            end
            
        else
            if ~(exist(energy.file{k}, 'file') == 2)
                error('Energy file does not exist')
            end
            load(energy.file{k},'En')
            EE = sort(En.Energies);
            Groundstate(k,:) = EE(1:n_groundstates);
            disp(['Vg: ', num2str(Vg_vec(k)), ', Groundstate: ', num2str(min(En.Energies)), ', first energy cut: ', ...
                num2str(sum(EE < (EE(1) + energy.first_cut))), ' states, second energy cut: ', num2str(sum(EE < (EE(1) + energy.second_cut))), ' states'])
            energy.energy_cut(k) = Groundstate(k,1) + energy.first_cut;
            energy.energy_cut_2(k) = Groundstate(k,1) + energy.second_cut;
            
        end
        
        
        
    end
    % plot(Vg_vec, Groundstate)
    
end

% ------------------------------ save -------------------------------------

parameters = vars2struct('condor', 'Vg_vec', 'Vb_vec', 'geometry', 'basis',...
                         'grid_mesh', 'energy', 'system', 'contact', 'bath', 'Beta', 'Vg_ran_cell', 'Vb_ran_cell');

path_to_save = './data/';
save([path_to_save, 'setup_', setup_id ,'.mat'],'parameters')

%check compile date

if ~isdeployed()
    change_list = dir([GROUNDPATH, 'exec_bm.m']);
    compile_list = dir([GROUNDPATH, 'exec_bm']);
    if isempty(compile_list)
        msgbox('Warning: Compiled program not found - please compile file, running compile.m in the bm_program folder', 'Compiled?', 'warn')
    else
        if datetime(change_list.datenum,  'ConvertFrom','datenum') > datetime(compile_list.datenum, 'ConvertFrom','datenum')
        msgbox('Warning: There were changes after compiling', 'Compiled?', 'warn')
        end
    end
    disp(['Compile date: ', compile_list.date])
end
disp(['Created file name: ', path_to_save, 'setup_', setup_id ,'.mat', ', number of Condor runs: ', num2str(numel(Vb_ran_cell))])

