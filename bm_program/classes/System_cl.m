classdef System_cl
    %author: Gerhard Dorn
    %date: august 15
    %use: create fermionic Hamiltonian, and solve it to produce object Energy
    %features: - 
    %          - 
    %          - 
    %general description of the system
    %   a subclass shall incorporate special features like conserving
    %   particle number and spin
    
    properties
        Parameter_list % form of table, depending on model type
        Model_name  %name like hubbard, kitaev, ...
        Basis %Basis used to set up many body Hamiltonian
        Hamiltonian_mb % in many body Basis (000, 001, 010, 100, 011, 101, ...)
        Tab
        
        %Hamiltonian_sp %either one sparse matrix or with block structure (single particle / system site, spin)
    end
    
    methods
        %constructor, depending on model, parameter handling has to be
        %consideredds
        function obj = System_cl(parameters, model, Basis, Hamiltonian_mb, Tab)
            obj.Parameter_list = parameters;
            obj.Model_name = model;
            obj.Hamiltonian_mb = Hamiltonian_mb;
            obj.Basis = Basis;
            obj.Tab = Tab;
            %obj.Hamiltonian_sp = 0;
            
            

        end
        
        %solve Hamiltonian and create object Energy (save outside of system)
        function Energ = f_solve(obj, attributes, energy_spacing)
            if nargin < 2, attributes.tolerance = 10^-12;    end
            field_names = fieldnames(attributes);
            if ~ismember('tolerance', field_names), attributes.tolerance = 10^-15; end
            if ~ismember('max_number', field_names), attributes.max_number = 20; end
            if nargin < 3 || isempty(energy_spacing), energy_spacing = inf; end
            
            %NOTE: energy_spacing is the range from the groundstate up to which eigenvalues will be calculated
            %(energy_cut_2)
            
            %NOTE: for spinless systems spin and spin_var are set to empty set []
            
            %TODO: solve Hamiltonian in each particle/spin sector if those
            %quantities are preserved anyway - procedure to check conservation
            
            %TODO: implement an approximate scheme for the relevant sectors - lanczos
            % also implement a version to have Green's functions representation via band lanczos 
            % implicit calculation of Q matrices
            
            %----> put spin_particle_conserved to properties???
            
            if strcmp(obj.Model_name, 'hubbard')
                spin_particle_conserved = true;
            else
                spin_particle_conserved = false;
            end
            
            %particle / spin conserved approach
            % -------------------------------
            % precondition: basis sorted in such a way
            
            %TODO: doesn't work for Heisenberg system! Docc = 0; implement new subclass?
            %TODO: 
            
            if spin_particle_conserved %energies will be block_diagonal in particle and spin sectors / intrinsic due to eigenvector
                % also works for spinless systems
                
                index = 0;
                %eigen_vector = sparse(obj.Basis.Numel_Bas, obj.Basis.Numel_Bas);
                %Tab = Table(obj.Basis,'NS');
                numel_block = obj.Tab.Numel_Blocks;
                eigen_vec_block = cell(numel_block,1);
                block_index = 0;
                
                %NOTE: preinitialize with inf
                [eigen_ener, part_num] = deal(inf(obj.Basis.Numel_Bas,1));
                if obj.Basis.N_spin > 1
                    spin = inf(obj.Basis.Numel_Bas,1);
                    structure = 'NSE';
                else
                    spin = [];  % for spinless systems
                    structure = 'NE';
                end
               
                
                
                for N_el = 0:obj.Basis.N_site*obj.Basis.N_spin
                    
                    if obj.Basis.N_spin == 1
                        LL = obj.Basis.Index('N',N_el); 
                        
                        if ~isempty(LL)
                            index = index(end) + (1:numel(LL));
                            block_index = block_index + 1;
                            % NOTE:  if this block is too large, calculate only the lowest #attributes.max_number
                            % defines how many eigenvalues per block will be solved:
                            if floor(numel(LL)/10)*10 > attributes.max_number
                                
                                % NOTE: Energy_spacing checks if we have calculated all eigenvalues up to the energy_cut_2 
                                % NOTE: Energy_spacing overrules max_number 
                                if energy_spacing < inf
                                    
                                    [eigen_vec_temp,eigenval] = eig(full(obj.Hamiltonian_mb(LL,LL)));
                                    %NOTE: not take the max_number lowest, but diagonalize the Hamiltonian full and
                                    %take as many needed to fullfil the energy_spacing criteria (not store all even
                                    %if full calculation can be done -> saves storage!!!
                                    %TODO: optimize, make recursivly
                                    L = diag(eigenval) < energy_spacing + min(diag(eigenval));
                                    eigen_vec_temp = eigen_vec_temp(:,L);
                                    eigenval = eigenval(L,L);
                                    disp(['N_el: ', num2str(N_el), ...
                                        ', numel_eigenvalues: ', num2str(sum(L)) ' of ' num2str(length(eigen_vec_temp))])
                                else
                                    k = attributes.max_number;
                                    [eigen_vec_temp, eigenval] = eigs(obj.Hamiltonian_mb(LL,LL),k, 'sa');
                                end
                                eig_length = size(eigen_vec_temp,1);
                                numel_eigen_vec = size(eigen_vec_temp,2);
                                % refill block with sparse zeros
                                eigen_vec = [sparse(eigen_vec_temp), sparse(eig_length,eig_length-numel_eigen_vec)];
                            else
                                [eigen_vec, eigenval] = eig(full(obj.Hamiltonian_mb(LL,LL)));
                            end
                            %NOTE: sets the first values ov eigen_ener in the current block to the calculated
                            %values
                            eigen_ener(index(1:numel(diag(eigenval)))) =  round(diag(eigenval), -floor(log10(attributes.tolerance)));
                            %NOTE: saves the eigenvector in a block (size(numel(LL), numel(LL))
                            eigen_vec_block{block_index} = eigen_vec;
                            
                            part_num(index) = N_el;
                        end
                    else
                        for spin_sec = obj.Basis.f_spin_range(N_el)
                            
                            %TODO: check if there are Basis sys
                            
                            LL = obj.Basis.Index('NS',N_el,spin_sec); %TODO: reevaluate - maybe shorter verion possible (with new functions)
                            
                            if ~isempty(LL)
                                index = index(end) + (1:numel(LL));
                                block_index = block_index + 1;
                                 % NOTE:  if this block is too large, calculate only the lowest #attributes.max_number
                                 % defines how many eigenvalues per block will be solved:
                                if floor(numel(LL)/10)*10 > attributes.max_number
                                    
                                    % NOTE: Energy_spacing checks if we have calculated all eigenvalues up to the energy_cut_2 
                                    % NOTE: Energy_spacing overrules max_number 
                                    if energy_spacing < inf
                                        
                                        [eigen_vec_temp,eigenval] = eig(full(obj.Hamiltonian_mb(LL,LL)));
                                        %NOTE: not take the max_number lowest, but diagonalize the Hamiltonian full and
                                        %take as many needed to fullfil the energy_spacing criteria (not store all even
                                        %if full calculation can be done -> saves storage!!!
                                        %TODO: optimize, make recursivly
                                        L = diag(eigenval) < energy_spacing + min(diag(eigenval));
                                        eigen_vec_temp = eigen_vec_temp(:,L);
                                        eigenval = eigenval(L,L);
                                        disp(['N_el: ', num2str(N_el), ', spin: ', num2str(spin_sec), ...
                                            ', numel_eigenvalues: ', num2str(sum(L)) ' of ' num2str(length(eigen_vec_temp))])
%                                         for k = [attributes.max_number:attributes.max_number:numel(LL), numel(LL)]
%                                             a_timer = tic;
%                                             eigenval = eigs(obj.Hamiltonian_mb(LL,LL), k, 'sa');
%                                             eig_time = toc(a_timer);
%                                             disp([num2str(k), ' eigenvalues: ', num2str(eig_time)])
%                                             if max(eigenval) > energy_cut
%                                                 disp(['Found maximal eigenvalue: ', num2str(max(eigenval)), ', number of eigenvalues necessary: ', num2str(k)])
%                                                 break
%                                             end
%                                             
%                                             if k > 200,
%                                                 disp('calculate full')
%                                                 
%                                             end
%                                             
%                                         end
                                    else
                                        k = attributes.max_number;
                                        [eigen_vec_temp, eigenval] = eigs(obj.Hamiltonian_mb(LL,LL),k, 'sa');
                                    end
                                    eig_length = size(eigen_vec_temp,1);
                                    numel_eigen_vec = size(eigen_vec_temp,2);
                                    % refill block with sparse zeros
                                    eigen_vec = [sparse(eigen_vec_temp), sparse(eig_length,eig_length-numel_eigen_vec)];
                                else
                                    [eigen_vec, eigenval] = eig(full(obj.Hamiltonian_mb(LL,LL)));
                                end
                                
                                eigen_ener(index(1:numel(diag(eigenval)))) =  round(diag(eigenval), -floor(log10(attributes.tolerance)));
                                
                                eigen_vec_block{block_index} = eigen_vec;
                                
                                spin(index) = spin_sec;
                                part_num(index) = N_el;
                            end
                            
                        end
                    end
                end

                % TODO: store eigen_vector also in block -> change settings in energy class for this 
                eigen_vector = obj.Tab.f_matrix(eigen_vec_block);
                
                [ener_unique,~,ind_position] = unique(eigen_ener);
                frequency = histc(eigen_ener(:),ener_unique);
                degeneracy = frequency(ind_position);
                if obj.Basis.N_spin > 1
                    spin_var = zeros(size(spin));
                else
                    spin_var = []; % for spinless systems
                end
                part_var = zeros(size(part_num));
                
                %------------------------------
                % general approach:
            else
                
                % full solve (develop also for groundstate only - eigs)
                [eigen_vec, eigenval] = eig(full(obj.Hamiltonian_mb));
                
                eigen_ener = round(diag(eigenval), -floor(log10(attributes.tolerance)));
                
                [ener_unique,~,ind_position] = unique(eigen_ener);
                frequency = histc(eigen_ener(:),ener_unique);
                degeneracy = frequency(ind_position);
                
                % -------determine particle sector and particle conservation ---------
                %determine if it is spin/particle conserved by evaluating the variance 
                %calculate the variance for particles (add + 1 to be able to use zero as not counting); var(x) = <x^2> - <x>^2
                Particle_distribution = (1 + repmat(obj.Basis.N_part,1, obj.Basis.Numel_Bas)).*(abs(eigen_vec) > attributes.tolerance);
                part_var = ( sum(Particle_distribution.^2,1)./sum(Particle_distribution~=0) - (sum(Particle_distribution,1)./sum(Particle_distribution ~= 0)).^2 ).';
                
                part_num = (sum(repmat(obj.Basis.N_part,1, obj.Basis.Numel_Bas).*(abs(eigen_vec) > attributes.tolerance),1)./ sum(abs(eigen_vec) > attributes.tolerance,1))';
                

                % -------determine spin sector and spin conservation ---------(if spin is defined)
                %calculate the variance for Spins (add + 0.5 to be able to use zero as not counting); var(x) = <x^2> - <x>^2
                if obj.Basis.N_spin > 1
                    Spin_distribution = (0.5 + repmat(obj.Basis.Spin,1, obj.Basis.Numel_Bas)).*(abs(eigen_vec) > attributes.tolerance);
                    spin_var = ( sum(Spin_distribution.^2,1)./sum(Spin_distribution~=0) - (sum(Spin_distribution,1)./sum(Spin_distribution~=0)).^2 ).';
                    
                    
                    spin = (sum(repmat(obj.Basis.Spin,1, obj.Basis.Numel_Bas).*(abs(eigen_vec) > tolerance),1)./ sum(abs(eigen_vec) > attributes.tolerance,1))';
                end
                
               
                %[eigen_ener, part_num, spin, degeneracy]
                
                %question of how to sort Energies???
                %analyse according to symmetries???
                
                %-------------------------------------
                % sort according to particle, spin, energies 
                % TODO: maybe change later to adopt different Hamiltonians, depends on situations
                %TODO: distinguish between particle, particle spin block structure (via variance)
                
                if obj.Basis.N_spin > 1
                    if var(part_var) < attributes.tolerance
                        if var(spin_var) < attributes.tolerance
                            
                            [~, ix] = sortrows([part_num, spin, eigen_ener], [1,2,3]);
                            structure = 'NSE';
                        else
                            [~, ix] = sortrows([part_num, eigen_ener], [1,2]);
                            structure = 'NE';
                        end
                    elseif var(spin_var) < attributes.tolerance && part_var > attributes.tolerance
                        [~, ix] = sortrows([spin, eigen_ener], [1,2]);
                        structure = 'SE';
                    else
                        ix = 1:numel(eigen_ener);
                        structure = 'E';
                    end
                    spin = spin(ix);
                else
                    if var(part_var) < attributes.tolerance
                        [~, ix] = sortrows([part_num, eigen_ener], [1,2]);
                        structure = 'NE';
                    else
                        ix = 1:numel(eigen_ener);
                        structure = 'E';
                    end
                    %NOTE: signalizes spinless systems
                    spin = [];
                    spin_var = [];
                end
                part_num = part_num(ix);
                
                eigen_ener = eigen_ener(ix);
                eigen_vector = eigen_vec(:,ix);
                degeneracy = degeneracy(ix);
            end

            
            Energ = Energy(eigen_ener, eigen_vector, part_num, part_var, spin, spin_var, degeneracy, obj.Basis, structure);
        end
        
        %visualizes system: either parameters, plot of Hamiltonian, spectrum,
        %some special characteristics
        function [] = s_print(obj)
            system_title = [obj.Model_name, ' model'];
            
            
            plot_settings;
            pos = get(gcf,'position');
            sizetitle = 15; sizelabel = 13;
            % headline
            annotation('textbox',  'units','pixel','position', [10, pos(4)*0.9, pos(3)-20,30], 'string', system_title, 'EdgeColor', 'none', ...
                'HorizontalAlignment', 'center','fontsize', sizelabel , 'interpreter', 'latex')
            
            %disp Hamiltonian
            [total_hamiltonian, parameters] = obj.s_hamiltonian;
            annotation('textbox',  'units','pixel','position', [10, pos(4)*0.7, pos(3)-20,30], 'interpreter', 'latex', 'string', total_hamiltonian, 'EdgeColor', 'none', ...
                'HorizontalAlignment', 'center','fontsize', sizelabel )
            
            %disp parameters
            
            annotation('textbox',  'units','pixel','position', [10, pos(4)*0.5, pos(3)/2-20,30], 'interpreter', 'latex', 'string', parameters(:,1), 'EdgeColor', 'none', ...
                'HorizontalAlignment', 'center','fontsize', sizelabel )
            
            if size(parameters,2) == 2
                 annotation('textbox',  'units','pixel','position', [10+pos(3)/2, pos(4)*0.5, pos(3)/2-20,30], 'interpreter', 'latex', 'string', parameters(:,2), 'EdgeColor', 'none', ...
                'HorizontalAlignment', 'center','fontsize', sizelabel )
            end
            
            %TODO: display matrix parameters or even in combination with plot of Geometry
            
            
            
        end
        
        %function to get a signature of the system for plots, filenames
        %etc...
        function info = s_disp(obj)
            
            
        end
        
        
        function [total_hamiltonian, parameters] = s_hamiltonian(obj)
            total_hamiltonian = '';
            parameters = {''};
        end
        
        %returns symmetries, conserved quantities, maybe analyses low
        %energy sector to implement 'low energy projection procedure' to
        %gain an effective Hamiltonian
        % -> analyse conservation of particle and spin sector already done in Solve
        function [] = Analyse(obj)
            
            
            
        end
        
        %maybe a function to adjust single parameter if influence is small
        function [] = change_System(parameters)
        
            
        end
    end
    
end

