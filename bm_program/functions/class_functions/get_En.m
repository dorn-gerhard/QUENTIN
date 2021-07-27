function En = get_En(setup_file_name, Vg)

if nargin < 2 || isempty(Vg), Vg = 0; end

load(setup_file_name, 'parameters');
struct2vars(parameters)






Bas = Basis(basis.N_site,basis.N_sector,basis.S_sector,basis.Docc_input,basis.N_spin); % class
Geom = Geometry(basis.N_site, geometry.structure, ones(basis.N_site)); % class



Hub = Hubbard(Geom, system.b + diag(contact.gate_coupling*Vg), system.U, system.xi, system.V, Bas, system.part_hole_sym); % class


En = Hub.f_solve(energy.calc, energy.second_cut); % class
%Q-matrix preinitialize


%En = En.f_Qmatrix_generate( energy.energy_cut_2(k));


