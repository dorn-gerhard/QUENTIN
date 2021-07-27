function master_equation(system, bath, contact, approximation)


%Note: first make a clear redefinition of all classes
% Basis, System_cl (Hubbard, ...), contact, Bath, Energy
% define how to hand over classes (filename if stored or nice parameter
% setup)

% define batch / slurm parameter research - discuss issue with eval -
% alternatively make a setup file for each run and delete everything in the
% end

%create constructors for simpler systems to ease the setup definition

% make a map and a design document for everything!!!
energy.basis
energy.system
energy.class = 'Energy.sf_from_system';

%Bas only for Hubbard needed and N_spin as property of Energy

%GF_Grid only for Bath needed and for CPT stuff



En = Energy.
system.class = 'Energy.sf_from_system';
system.basis = basis;
system.hamiltonian = ; %previous system
En = get_class(system);

En = Energy.sf_from_system(system

%think about how to use Energy Threshold 1 and Energy Threshold 2!!!


  if ~isempty(approximation.partial_secular_approximation)
        En = condense_energy_blocks(En, approximation.partial_secular_approximation);
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

%define energy_cut 1 and energy_cut 2 from approximation and En!!

if En.Basis.N_spin == 2 && isempty(energy.spin_tag)
       %TODO: check out if energy is used!
   En = En.f_resort('NE'); 
   %(for Basis creation - make unique if spin is used or not!!! - but
   %where? - in System_cl.f_solve (currently only implemented for conserved
   %spin and particle number
end

        [K, I, Sig, ~, Tab_cut, status, calc_time, column_sum, ev_min, choi_value, X_cur, H_LS_calculated] = lindblad(En, ...
            contact.contact_site, Tab, Mu, contact.coupling, Beta, bat, condor, Energy_cut, Energy_cut_2, run_id, min_en_ind, energy.eigenvalue_tolerance, multi_symmetrization);
