%test_bath file for all possible combinations


GROUNDPATH = '/temp/dorn/run/bm_program/';

	addpath(GROUNDPATH)
	addpath([GROUNDPATH, 'classes/'])
	addpath([GROUNDPATH, 'functions/'])
	addpath([GROUNDPATH, 'functions/class_functions/'])
    addpath([GROUNDPATH, 'functions/c_functions/m_GF/'])
	addpath([GROUNDPATH, 'functions/c_functions/'])
	addpath([GROUNDPATH, 'functions/plot/'])
    addpath([GROUNDPATH, 'functions/test_functions/'])



%a) spin conserved, spin=2, spin symmetric, 2 +1 orbitals

%b) spin_conserved, spin = 2, not symmetric 3 orbitals (different temperatures for different spin)

%c) not spin conserved, spin = 2, mixing of up and down spin by angle (polarization)

%d) spinless, 3 baths, 2+1+3 orbitals

type = 'a';


%basis

basis.N_site = 4;
basis.N_sector = []; %empty for all
basis.S_sector = []; %empty for all
basis.Docc_input = 1; % 1 for all, 0 for up and down, 2 for none, up and down; 3 for up, down and double

%grid
grid_mesh.uniform = false;
grid_mesh.lim = [-6,6];
grid_mesh.numel_points = 24001; %24001 
grid_mesh.imag_i0plus = 10^-7;
grid_mesh.points = [-20:-10, -9.5:0.5:-6.5, linspace(-6,6,19345), 6.5:0.5:9.5, 10:20];

if  grid_mesh.uniform == false
    GF_grid = Grid(grid_mesh.points, grid_mesh.imag_i0plus);
else
    GF_grid = Grid(linspace(grid_mesh.lim(1), grid_mesh.lim(2) ,grid_mesh.numel_points), grid_mesh.imag_i0plus); % class
end

%bath

bath.type = 'tight_binding';%'wide_band'; %'from_file'; % 
bath.integrator = 'Fast_Detail'; %'Fast', 'Matlab', 'Detail', 'Constant'

% for wide band and tight binding leads

if type == 'a'
    %a) spin conserved, spin=2, spin symmetric, 2 +1 orbitals
    
    bath.N_leads = 2;
    basis.N_spin = 2;
    bath.spin_conserved = [true; true]; %mixture of spins (mostly true) - make an example for spin_conserved = false
    bath.spin_symmetric = [true; true]; %same bath for spin up and down (spin_conserved necessary)
    bath.numel_orbitals = [2   ; 1]; %should correspond to t_leads and eps_0 (if it is half then spin is not conserved)
    bath.shift_bath =     [true; false];
    
    bath.t_leads =  {[5, 1; 1, 3]; 8};
    bath.eps_0 = {[-1, 0; 0, 2]; 3};
    beta = [30;20];
    mu = [-2;2];
    
    contact.coupling  = {[1,1; 0,0; 0,0; 0,0];...
        [  0;   0;   1;   1]};
    
elseif type == 'b'
    %b) spin_conserved, spin = 2, not symmetric 3 orbitals (different temperatures for different spin)
    
    bath.N_leads = 3;
    basis.N_spin = 2;
    bath.spin_conserved = [true; true; true]; %mixture of spins (mostly true) - make an example for spin_conserved = false
    bath.spin_symmetric = [false; false; false]; %same bath for spin up and down (spin_conserved necessary)
    bath.numel_orbitals = [2   ; 1; 1]; %should correspond to t_leads and eps_0 (if it is half then spin is not conserved)
    bath.shift_bath =     [true; false; false];
    
    bath.t_leads =  {[5, 1 0 0; 1, 3, 0, 0; 0, 0, 2, 1; 0, 0, 1, 2];...
                      [8 0; 0 0]; ...
                      [0 0; 0 2]};
    bath.eps_0 = {[-1, 0, 0, 0; 0, 2, 0, 0; 0, 0, 0, 0; 0, 0, 0, 1];...
                   [3, 0; 0, 0]; ...
                   [0, 0; 0, 1]};
    beta = [30;20; 40];
    mu = [-2; 2; 2];
    
    contact.coupling  = {[1,1,0,1; 0,0,0,0; 0,0,0,0; 0,0,0,0];...
        [  0, 0;   0,0;   1,1;   1,1]; ...
        [  0, 0;   0,0;   0.8, 0.8; 1,1]};
 
    
    
    
elseif type == 'c'
    %c) not spin conserved, spin = 2, mixing of up and down spin by angle (polarization)

    bath.N_leads = 2;
    basis.N_spin = 2;
    bath.spin_conserved = [false; true]; %mixture of spins (mostly true) - make an example for spin_conserved = false
    bath.spin_symmetric = [false; true]; %same bath for spin up and down (spin_conserved necessary)
    bath.numel_orbitals = [2   ; 2]; %should correspond to t_leads and eps_0 (if it is half then spin is not conserved)
    bath.shift_bath =     [true; false];
    
    bath.t_leads =  {[5, 1 3 2;...
                     1, 3, 1, 0; ...
                     3, 1, 2, 1; ...
                     2, 0, 1, 2];...
                                     [8 1; ...
                                     1 2]};
    bath.eps_0 = {[-1, 0, 0, 0;...
                    0, 2, 0, 0;...
                    0, 0, 0, 0;...
                    0, 0, 0, 1];...
                                    [3, 0; ...
                                    0, 0]};
    beta = [30;20];
    mu = [-2; 2];
    
    contact.coupling  = {[1,1,0,1; ...
                          0,0,0,0; ...
                          0,0,0,0; ...
                          0,0,0,0];...
                                    [0, 0;...
                                     0, 0;...
                                     1, 1; ...
                                     1, 1]};
 
    
elseif type == 'd'
    %d) spinless, 3 baths, 2+1+3 orbitals

    
end

% bath.t_leads = cell(bath.N_leads, 1); %rows: baths, columns: spin
% bath.eps_0   = cell(bath.N_leads,  1); %rows: baths, columns: spin
% [bath.t_leads{:}] = deal(t_leads); %same bath on both sides
% [bath.eps_0{:}] = deal(eps_0);



bath_temp = Bath.c_tight_binding(GF_grid, bath.t_leads, bath.shift_bath, bath.eps_0, bath.integrator, basis.N_spin, bath.numel_orbitals, beta, mu);
%plot_GF(bath_temp.GF{1}, GF_grid)
%coupling

contact.couple_bias_gate = false;
contact.contact_site = [1:basis.N_site];
contact.gate_coupling = ones(1,basis.N_site);



Bath.f_check_Coupling(bath_temp, contact.coupling, basis)


bath_temp.Mu(1) = -5

bath_temp.s_plot_DOS(contact.coupling)
