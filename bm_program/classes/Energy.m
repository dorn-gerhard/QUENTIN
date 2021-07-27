classdef Energy < handle
    %author: Gerhard Dorn
    %date: august 15
    %use: create fermionic Hamiltonian, and solve it to produce object Energy
    %features: -
    %          -
    %          -
    %general description of the isolated central system (no influence of bath)
    %   a subclass shall incorporate special features like conserving
    %   particle number and spin
    
    properties
        %Energies        % full solution of the System
        
        %Degeneracy      % frequency of each Eigenenergy
        Basis           % Basis in which the Eigenstates are expressed
        %Index           % Index
        Qmatrix         % Qmatrix system in Particle, Spin block, all stored, use g_Qmat to load wished Q_matrix, works also without spin
        Qmatrix_counter % counts number of Q matrices
        
        Q_table         % Qtable corresponding to Energies (old)
        Eigenvector     % corresponding Eigenstates of the System
        % special Properties of the solution of a Hamiltonian that preserves Spin and Particle Number
        % thus eigenspaces have distinct spin and part. numb. -> should be part of a subclass
        
        %N_part          % particle number of Eigenstate
        %Spin            % Spin number of Eigenstate
        %N_part_var      % particle variance (0 for all if particle number is conserved)
        %Spin_var        % Spin variance (0 for all if spin is conserved)
        Structure       % Structure, represents in which order/grouping eigenstates are arranged ('NSE')
        Numel_En        % Number of eigenstates
        Categ
        DBUG = true;
    end
    
    
    methods
        %list:
        % Energy ...                constructor
        % f_c ...                   creates Q_matrix in special particle and spin sector
        % f_Qmatrix_generate ...    creates all Q_matrices for all sectors and all possible system sites and spins
        % g_Qmat ...                gets Qmatrix of corresponding particle and spin sector according to site, spin and dagger
        % f_create_Q_table ...      creates Qtable class type object (OLD - use Q_Table class instead)
        % s_site_info ...           shows conversion to system sites
        % f_En ...                  function to return Energies, Eigenvectors - very dynamical
        % g_En ...                  getter very plain to return Energies
        % f_Ev ...                  function to return Eigenvecotrs - very dynamical
        % f_In ...                  function to return index
        % s_print ...               shows a table of energies
        
        %constructor - somehow incorporate System
        function obj = Energy(eigen_ener, eigen_vec, part_num, part_var, spin, spin_var, degeneracy, basis, structure)
            % obj = Energy(eigen_ener, eigen_vec, part_num, part_var, spin, spin_var, degeneracy, basis, structure)
            %
            % Energy  - Constructor
            %   Energy(eigen_ener, eigen_vec, part_num, spin, degeneracy, basis)
            %   creates Energy object as full solution of the Hamiltonian (System class)
            %   Eigenvalues / vectors are stored in particle / Spin and Energy increasing order
            %   in the basis representation stored in the basis class
            %
            obj.Numel_En = numel(eigen_ener);
            
            if isempty(spin_var) || numel(spin_var) ~= obj.Numel_En, spin_var = []; end
            if isempty(part_var) || numel(part_var) ~= obj.Numel_En, part_var = []; end
            
            %obj.Energies = eigen_ener;
            %TODO: store Eigenvector in Block, corresponding to structure - reprogramm the following functions
            % --- g_EV
            % --- f_EV (old - to be deleted)
            % --- f_c_elem (transfere to Eigen_Vector)
            % --- f_site_info
            % --- Eigen_Vector
            % maybe keep g_Ev and Eigen_Vector (one for fast search, one for detailed search
            % make Eigenvector private property
            
            obj.Eigenvector = eigen_vec;
            
            
            %obj.N_part = part_num;
            %obj.Spin = spin;
            %obj.N_part_var = part_var;
            %obj.Spin_var = spin_var;
            %obj.Degeneracy = degeneracy;
            obj.Basis = basis;
            
            %obj.Index = (1:numel(obj.Energies)).';

            %obj.Q_table = obj.f_create_Q_table();   %OLD, TODO: still used in SCB_Y_mat, SCB_algo and SCB_test_integral, SCB_Y_mat_new
            %obj.Structure = structure;
            
            % TODO: maybe reevaluate storage of Energies, Eigenvector, Spin + var, N_part + var, Degeneracy, 
            % option to organize in Category (see below) -> manage access
            % option to include not full System Energies
            
            obj.Categ = Category({'Index', 'I'}, (1:obj.Numel_En).');
            
            obj.Categ = obj.Categ.f_append({'Energies', 'E'}, eigen_ener, 10^-10);
            %NOTE: Define Spin in Categ class only if there is a spin in the basis!
            if basis.N_spin > 1
                obj.Categ = obj.Categ.f_append({'Spin', 'S'}, spin, [], spin_var);
            end
            obj.Categ = obj.Categ.f_append({'N_part', 'N'}, part_num, [], part_var);
            obj.Categ = obj.Categ.f_append({'Degeneracy', 'D'}, degeneracy);
            %NOTE: N_part, Spin (for spin systems) give average numbers if part_var and spin_var are not zero!!!
            % the block structure feature is stored in structure
            
            obj.Structure = structure; %obj.Categ.f_shortnames();

            %obj.Qmatrix{dagger/2+1.5, N_part_inp + 1, Spin_index, c_site, c_spin_index}
            
            %TODO: distinguish between spinless and spin systems!!!
            
            %NOTE: Qmatrix indices:
            %1) no dagger, dagger [2]
            %2) particle sector [basis.N_spin*basis.N_site+1] (if particle number is conserved)
            %3) spin sector [2*N_site + 1] (if spin is conserved (indepedendent of conserved particle number)
            %4) acting site of creation / annihilation operator
            %5) acting spin of -""- [2] (if not spinless)
            if obj.Basis.N_spin == 1
                spin_sector_range = 1;
                operator_spin_range = 1;
                if ismember('N', obj.Structure)
                    particle_sector_range = basis.N_site + 1;
                else
                    particle_sector_range = 1;
                end
                
            else
                operator_spin_range = 2;
                if ismember('N', obj.Structure)
                    particle_sector_range = 2*basis.N_site + 1;
                else
                    particle_sector_range = 1;
                end
                if ismember('S', obj.Structure)
                    spin_sector_range = 2*basis.N_site + 1;
                else
                    spin_sector_range = 1;
                end
            end
            
            obj.Qmatrix = cell(2, particle_sector_range, spin_sector_range, basis.N_site, operator_spin_range);
            obj.Qmatrix_counter = 0;
            
            %{
            % TODO: check if this makes sense - Eigenenergy and Vector will always be defined
            % N_part and Spin + variances not always...
            
            obj.Category(sum_indices  + 2).Name = 'Vector'; 
            obj.Category(sum_indices  + 2).Value = eigen_vec;
            obj.Category(sum_indices  + 2).Tol = 0;
            
            obj.FILTER_CAT_NAMES = {obj.Category(:).Name}; %full name
            obj.Category(1).Structure = cellfun(@(x) upper(x(1)), obj.FILTER_CAT_NAMES);
            
            keyset = {'Index', 'Energies', 'N_part', 'Spin', 'Degeneracy', 'Eigenvectors', 'Squared Spin'};
            valueset = {1, E_index + 1, N_index + 1 , S_index + 1, sum_indices + 1, sum_indices + 2, sum_indices + 3};
            obj.FILTER_CAT_NAMES = containers.Map(keyset, valueset);
            %TODO: make a list which contains all category names, so that it can be used for referencing
            %}
        end
        
        
        function index = Index(obj, varargin)
            index =  obj.Categ.f_access('Index', varargin{:});
        end
        
        function index = Gap_Energies(obj, varargin)
            index =  obj.Categ.f_access('Gap_Energies', varargin{:});
        end
        
        function energ = Energies(obj, varargin)
            energ =  obj.Categ.f_access('Energies', varargin{:});
        end
        
        function n_part = N_part(obj, varargin)
            n_part =  obj.Categ.f_access('N_part', varargin{:});
        end
        
        function spin = Spin(obj, varargin)
            spin = obj.Categ.f_access('Spin', varargin{:});
        end
        
        function degeneracy = Degeneracy(obj, varargin)
            degeneracy = obj.Categ.f_access('Degeneracy', varargin{:});
        end
        
        function aver_ener = Average_Energies(obj, varargin)
            aver_ener = obj.Categ.f_access('Average_Energies', varargin{:});
        end
        
        function eigen_vector = Eigen_Vector(obj,varargin)
            %NOTE: Return just the eigenvectors of varargin (columns) but all relevant rows (dependent on NS
            %preservation)
            if nargin < 2 || isempty(varargin) 
                eigen_vector = obj.Eigenvector;
            else
                filter = obj.Categ.filter_inp(varargin{:});
                if ismember('E', varargin{1}) || ismember('A', varargin{1})
                    %NOTE: to prevent search in Basis.Categ with filter E or A which cannot be defined!
                    filter_names = varargin{1};
                    L = filter_names=='E' | filter_names == 'A';

                    filter2 = obj.Basis.Index(filter_names(~L), varargin{find(~L)+1});
                else
                    filter2 = obj.Basis.Index(varargin{:});
                end
                
                eigen_vector = obj.Eigenvector(filter2, filter);
            end
            
        end
        

        
        %NOTE: rewritten
        function obj = f_resort(obj, new_structure)
            
           
            % first find permutation 
            A = zeros(obj.Numel_En,0);
            ind = zeros(size(new_structure));
            for k = 1:numel(new_structure)
                ind(k) = obj.Categ.List(new_structure(k));
                A(:,k) = obj.Categ.Data(ind(k)).Values;
            end
           
           [~, index] = sortrows(A);
            
            % second apply permutation to Categ
            for k = 1:obj.Categ.Numel_Data_Types   
                % maybe exclude index!!!
                if obj.Categ.Data(k).Shortname ~= 'I'
                    obj.Categ.Data(k).Values = obj.Categ.Data(k).Values(index);
                end
            end
            
            
            % third apply permutation to Eigenvec
            % NOTE: columns belong to eigenvalues, rows to basis
            obj.Eigenvector = obj.Eigenvector(:,index);
            
            
            % fourth reset Qmatrix
            obj = obj.f_reset_Qmatrix;
            obj.Structure = new_structure;
            
        end
        
        
        %TODO: Outdated - maybe an unsort function is not wise
        function obj = f_unsort(obj)
            % obj = f_unsort(obj)
            % 
            % used to unsort previously sorted Eigenbasis according to Structure (hasn't been changed when sorted)
            [~, index] = sort(obj.Index);
            %TODO: change of Data has to be implemented in Category!!!
            obj.Energies = obj.Energies(index);
            obj.N_part = obj.N_part(index); 
            obj.Spin = obj.Spin(index);
            obj.Eigenvector = obj.Eigenvector(:,index);
            obj.Index = obj.Index(index,:);
            obj.N_part_var = obj.N_part_var(index);
            obj.Spin_var = obj.Spin_var(index);
            obj.Degeneracy = obj.Degeneracy(index);
            
            
            % work also in Categories (Energies, N_part, Spin, Index, Degeneracy)
            for n_Cat = 1: numel(Category)
                obj.Category(n_Cat).Values = obj.Category(n_Cat).Values(index);
            end
            
        end
        
        % apply operator - or implement directly as Q-matrix system
        % TODO: shall depend on structure (if spin is not conserved - ignore spin block - depending on structure)
        function [Q] = f_c(obj, dagger, N_part_inp, Spin_inp, c_site, c_spin, energy_cut)
            if nargin < 7, energy_cut = []; end
            % c  - method - returns Q matrix in eigenbasis from certain sectors (if the quantum numbers (particle /
            %   spin ) are conserved
            %   c(dagger, N_part_inp, Spin_inp, c_site, c_spin)
            %   N_part_inp and Spin_inp define the spin/particle sector of the ket (right) vectorspace
            %   c_{c_site,c_spin}^dagger is acting on
            %   Therefore bra (left) vectorspace is <N_part_inp + dagger, Spin_inp + dagger*c_spin|
            
            %NOTE: steps to calculate Q-matrix:
            % first:  go to local space basis and look what creation / annihilation operators are doing, create
            %         Matrix K which includes the Fermi sign
            % second: Multiply the corresponding eigenstates (eigenvectors) onto the K Matrix, to get the Q-Matrix
            %         in the eigenspace basis. Incorporate energy cuts.
            
            
            %NOTE: should distinguish between the following cases
            
            %1) spinless
            %1a) spinless + particle number conserved
            %               Spin_inp, c_spin = []
            %1b) spinless + particle number not conserved
            %               N_part_inp, Spin_inp, c_spin = []
            %
            %2) with spin
            %2a) spin + particle number and spin conserved
            %
            %2b) spin + particle number conserved
            %               Spin_inp = []
            %2c) spin + spin number conserved
            %               N_part_inp = []
            %2d) spin
            %               Spin_inp, N_part_inp = []
            
            %      N_part_inp  |  Spin_inp  |  c_spin
            % ___________________________________________
            % 1a)|     #       |     [ ]    |    [ ]
            % 1b)|    [ ]      |     [ ]    |    [ ]
            % 2a)|     #       |      #     |     #
            % 2b)|     #       |     [ ]    |     #
            % 2c)|    [ ]      |      #     |     #
            % 2d)|    [ ]      |     [ ]    |     #
            
            % Define modes
            if obj.Basis.N_spin == 1
                if ~isempty(Spin_inp) || ~isempty(c_spin)
                    warning('In spinless case, Spin_inp and c_spin should not be defined but empty!')
                end
                if ~ismember('N', obj.Structure)
                    mode = 2;
                    filter_arguments = {};
                    filter_names = '';
                    if ~isempty(N_part_inp)
                        warning('Case 1b, particle number not conserved, N_part_inp should not be defined but empty!')
                    end
                else
                    mode = 1;
                    filter_arguments = {N_part_inp};
                    filter_names = 'N';
                    if isempty(N_part_inp)
                        warning('Use g_full_Qmat or define N_part_inp')
                    end
                end
            else
                % with spin
                
                if ismember('S', obj.Structure)
                    %checks if right particle/spin sector is well defined
                    if all(~ismember(Spin_inp,obj.Basis.f_spin_range(N_part_inp)))
                        warning('right side: spin out of range')
                        disp(['N_right: ', num2str(N_part_inp), ', S_right: ', num2str(Spin_inp), ', Dag: ', num2str(dagger), ', c_spin: ', num2str(c_spin)])
                    end
                    if all(~ismember(Spin_inp + dagger*c_spin, obj.Basis.f_spin_range(N_part_inp + dagger)))
                        warning('left side: spin out of range')
                        disp(['N_left: ', num2str(N_part_inp + dagger), ', S_left: ', num2str(Spin_inp + dagger*c_spin), ', Dag: ', num2str(dagger), ', c_spin: ', num2str(c_spin)])
                    end
                    
                    
                    if ismember('N', obj.Structure)
                        
                        mode = 3;
                        filter_names = 'NS';
                        filter_arguments = {N_part_inp, Spin_inp};
                    else
                        mode = 5;
                        filter_names = 'S';
                        filter_arguments = {Spin_inp};
                    end
                else
                    if ismember('N', obj.Structure)
                        mode = 4;
                        filter_names = 'N';
                        filter_arguments = {N_part_inp};
                    else
                        mode = 6;
                        filter_names = '';
                        filter_arguments = {};
                    end
                end
            end
            
            if ismember('A', obj.Structure)
                energy_tag = 'A';
            else
                energy_tag = 'E';
            end
            
            % check input: left and right sector
            
            
            
            N_site = obj.Basis.N_site;
            
            %bas_bin = obj.Basis.f_bin_(N_part_inp);
            %bas = obj.Basis.f_N(N_part_inp);
            
            bas_bin = obj.Basis.f_bin(filter_names, filter_arguments{:});                % particle_sector dependent
            bas = obj.Basis.Bas_Dez(filter_names, filter_arguments{:});                 % particle_sector dependent
            
            if obj.Basis.Smaller_Binary_Sector_Used_For(1) == 'd'
                spindown_shift = N_site;
            else
                spindown_shift = 0;
            end
            spinup_shift = N_site - spindown_shift;
            
            if c_spin == -1        % mode 3,4,5,6
                spin_shift = spindown_shift;
            elseif c_spin == 1     % mode 3,4,5,6
                spin_shift = spinup_shift;
            elseif isempty(c_spin) % mode 1,2
                spin_shift = 0;
            end
            
            %Create a logical mask which indicates which basis state has an electron at c_site that can be created/annihilated
            %this is the column index in the Q_matrix
            LL = find(bas_bin(:,c_site + spin_shift) == (dagger-1)*(-1/2));
            %Calculate the Fermi sign (how many states are before that position?
            %Remember Basis setup: 1 0 0 1 0 1: c_6^t c_4^t c_1^t |0 >
            %Fermi will be the entry in the Q_matrix in local space basis
            Fermi = (-1).^(sum(bas_bin(LL,1:(c_site + spin_shift-1)),2));
            
            %calculate the new states that come out by applying the creation / annihilation operator (for
            %each state it is exactly one new state (count from the right)
            new = bas(LL) + dagger * 2^(obj.Basis.N_spin*N_site - (c_site + spin_shift));
            
            %next step: find index in new particle sector
            %NOTE: basis has always an N Category, the question is, if the particle number is conserved,
            % which is stored in Energy.Structure!
            if ismember('N', obj.Structure)
                if ismember('S', obj.Structure)
                    % mode 3
                    left_bas_total = obj.Basis.Bas_Dez('NS', N_part_inp + dagger, Spin_inp + dagger*c_spin);
                else
                    % mode 1, 4
                    left_bas_total = obj.Basis.Bas_Dez('N', N_part_inp + dagger);            % particle_sector dependent
                end
            else
                if ismember('S', obj.Structure)
                    % mode 5
                    left_bas_total = obj.Basis.Bas_Dez('S', Spin_inp + dagger*c_spin);
                else
                    % mode 2, 6
                    left_bas_total = obj.Basis.Bas_Dez();                                   % particle_sector dependent
                end
            end
            %left_bas_index is the row index in the Q_matrix
            [~, left_bas_index] = ismember(new, left_bas_total);
            %left_bas = left_bas_total(left_bas_index);
            
            %creating the sparse Q-matrix in local space basis
            K = sparse(left_bas_index, LL, Fermi, numel(left_bas_total), numel(bas));
            %K = sparse (numel(left_bas_total), numel(bas));
            %K(sub2ind(size(K),left_bas_index,LL)) = Fermi;
            
            % Q = obj.g_Ev(N_part_inp + dagger, Spin_inp + dagger * c_spin)' * K * obj.g_Ev(N_part_inp, Spin_inp);
            
            if ~isempty(energy_cut)
                % Provide a practical
                %NOTE: Altough the Category class would ignore filter names which are not defined, it
                % does not represent the preserved quantities!, N and S could be average numbers for each
                % eigenstate! Therefore distinguish!
                
                %create a logical subindex vector to know where to store the values
                %TODO: Replace find function for index_b, index_a by Categ.Subindex function (implemented also as Energy.Subindex(Filter, Arguments, Subfilter, Arguments))
                if ismember('N', obj.Structure)
                    if ismember('S', obj.Structure)
                        % mode 3
                        index_b = obj.Subindex('NS', {N_part_inp, Spin_inp}, energy_tag, {{-inf, energy_cut}});
                        index_a = obj.Subindex('NS', {N_part_inp + dagger, Spin_inp + dagger*c_spin}, energy_tag, {{-inf, energy_cut}});
                        
                        %index_b = find(obj.Energies('NS', N_part_inp, Spin_inp) <= energy_cut);
                        %index_a = find(obj.Energies('NS', N_part_inp + dagger, Spin_inp + dagger*c_spin) <= energy_cut);
                        temp = obj.Eigen_Vector(['NS', energy_tag], N_part_inp + dagger, Spin_inp + dagger*c_spin, {-inf, energy_cut})' * ...
                            K * obj.Eigen_Vector(['NS', energy_tag], N_part_inp, Spin_inp, {-inf, energy_cut});
                        
                        
                    else
                        % mode 1, 4
                        index_b = obj.Subindex('N', {N_part_inp}, energy_tag, {{-inf, energy_cut}});
                        index_a = obj.Subindex('N', {N_part_inp + dagger}, energy_tag, {{-inf, energy_cut}});
                        
                        %index_b = find(obj.Energies('N', N_part_inp) <= energy_cut);
                        %index_a = find(obj.Energies('N', N_part_inp + dagger) <= energy_cut);
                        temp = obj.Eigen_Vector(['N', energy_tag], N_part_inp + dagger, {-inf, energy_cut})' * ...
                            K * obj.Eigen_Vector(['N', energy_tag], N_part_inp, {-inf, energy_cut});
                    end
                else
                    if ismember('S', obj.Structure)
                        % mode 5
                        index_b = obj.Subindex('S', {Spin_inp}, energy_tag, {{-inf, energy_cut}});
                        index_a = obj.Subindex('S', {Spin_inp + dagger*c_spin}, energy_tag, {{-inf, energy_cut}});
                        
                        %index_b = find(obj.Energies('S', Spin_inp) <= energy_cut);
                        %index_a = find(obj.Energies('S', Spin_inp + dagger*c_spin) <= energy_cut);
                        temp = obj.Eigen_Vector(['S', energy_tag], Spin_inp + dagger*c_spin, {-inf, energy_cut})' * ...
                            K * obj.Eigen_Vector(['S', energy_tag], Spin_inp, {-inf, energy_cut});
                        
                    else
                        % mode 2, 6
                        index_b = obj.Index(energy_tag, {-inf, energy_cut});
                        index_a = obj.Index(energy_tag, {-inf, energy_cut});
                        
                        %index_b = find(obj.Energies <= energy_cut);
                        %index_a = find(obj.Energies <= energy_cut);
                        temp = obj.Eigen_Vector(energy_tag, {-inf, energy_cut})' * K * obj.Eigen_Vector(energy_tag,{-inf, energy_cut});
                        
                    end
                end
                
                %NOTE: index_aa shall change in column direction (since Matlab addresses matrix elements first down then right)
                [index_bb,index_aa] = meshgrid(index_b, index_a);
                %store final sparse Q_matrix (keep dimensions of block)
                Q = sparse(index_aa(:), index_bb(:), temp(:), numel(left_bas_total), numel(bas));
                
                %store the final sparse Q-matrix
                %Q = sparse(numel(left_bas_total), numel(bas));
                %Q(index_a, index_b) = temp;
                size_temp = [' effective size: ',num2str(size(temp))];
            else
                if ismember('N', obj.Structure)
                    if ismember('S', obj.Structure)
                        % mode 3
                        Q = obj.Eigen_Vector('NS', N_part_inp + dagger, Spin_inp + dagger * c_spin)' * K * obj.Eigen_Vector('NS',N_part_inp, Spin_inp);
                        
                    else
                        % mode 1, 4
                        Q = obj.Eigen_Vector('N',N_part_inp + dagger)' * K * obj.Eigen_Vector('N',N_part_inp);
                        % Q = obj.g_Ev(N_part_inp + dagger,[])' * K * obj.g_Ev(N_part_inp,[])
                    end
                else
                    if ismember('S', obj.Structure)
                        % mode 5
                        Q = obj.Eigen_Vector('S', Spin_inp + dagger * c_spin)' * K * obj.Eigen_Vector('S', Spin_inp);
                        
                    else
                        % mode 2, 6
                        Q = obj.Eigen_Vector()' * K * obj.Eigen_Vector();
                    end
                end
                
                size_temp = ' full calculation!';
            end
            obj.Qmatrix_counter = obj.Qmatrix_counter + 1;
            if ~isdeployed()
                disp(['Qmatrix number ', num2str(obj.Qmatrix_counter), ' created! ', num2str(size(Q)), size_temp])
                disp(['dag: ', num2str(dagger), ', c_site: ', num2str(c_site), ', c_spin: ', num2str(c_spin), ', N_part: ', num2str(N_part_inp), ', Spin: ', num2str(Spin_inp)])
            end
            
            
            
        end
        
        
        %TODO: make a general function that allows to calculate < N | c^dag_(site, spin) | M > for two distinct
        %eigenstates N and M
        function c_op = f_c_elem(obj,dagger,c_site,c_spin,index_a, index_b)
            
            % index_a describes an eigenstate - find appropriate sector to operate in (if particle spin symmetric
            % and so on, depending on obj.Structure, 
            
            
              if obj.Basis.N_spin == 1
                %spinless fermions, ignore Spin_inp and c_spin
                %TODO: fuse with functions below
                
                
                %zerost: check if energy structure inherits 'N' (do eigenstates live in a certain particle sector?)
                if ~ismember('N', obj.Structure)
                    error('not yet implemented - general implementation needed - cost intensive');
                end
                
                % check if < index_a |c| index_b > can work:
                if ~((obj.N_part(index_a) == obj.N_part(index_b) + dagger) & (obj.Spin(index_a) == obj.Spin(index_b) + dagger * c_spin))
                    warning('< a | c^dagger_{c_site c_spin} | b > is zero'),
                end
                % tactics:
                % first take the binary basis sector in which the first state lives
                % then apply the operator (mind the fermi sign according to binary convention)
                % then multiply with eigenstate from left and eigenstate from right

                %first: find in which spin/particle sector eigenstate is situated
                N_part_b = obj.N_part(index_b);
                Spin_b = obj.Spin(index_b);
                
                %bas_bin = obj.Basis.f_bin_(N_part_b); %Spin argument left empty, since there is no spin
                %bas = obj.Basis.f_N(N_part_b);
                bas_bin = obj.Basis.f_bin('N', N_part_b);
                bas = obj.Basis.Bas_Dez('N',N_part_b);
                
                EV_b =  obj.Eigenvector(:,index_b);
                
                N_site = obj.Basis.N_site;
                
                % could be deleted, but simpler to set to 0 
                spin_shift = 0;
                
                LL = find(bas_bin(:,c_site + spin_shift) == (dagger-1)*(-1/2)); %searches at which states in space basis the operator c can act (without getting 0) 
                Fermi = (-1).^(sum(bas_bin(LL,1:(c_site + spin_shift-1)),2));
                new = bas(LL) + dagger * 2^(obj.Basis.N_spin*N_site - (c_site + spin_shift));
                
                
                %left_bas_total = obj.Basis.f_N(N_part_b + dagger);
                left_bas_total = obj.Basis.Bas_Dez('N',N_part_b + dagger);
                [~, left_bas_index] = ismember(new, left_bas_total);
                %left_bas = left_bas_total(left_bas_index);
                K = sparse (numel(left_bas_total), numel(bas));
                K(sub2ind(size(K),left_bas_index,LL)) = Fermi;
                
                Lb = obj.Basis.f_ind(N_part_b, Spin_b);
                La = obj.Basis.f_ind(N_part_b + dagger, Spin_b + dagger*c_spin);
                EVa = obj.Eigenvector(La, index_a)' ;
                EVb = obj.Eigenvector(Lb, index_b);
                c_op = EVa * K * EVb;
              end
            
            
        end
        
        %works for Spin and Spinless, used for local storage in cell (fast)
        function obj = f_Qmatrix_generate(obj, energy_cut)
            if nargin < 2 || isempty(energy_cut), energy_cut = []; end
            %NOTE: distinguish between 6 modes (spinless, with spin, and so on... )
            
            %NOTE: Qmatrix indices:
            %1) no dagger, dagger [2]
            %2) particle sector [basis.N_spin*basis.N_site+1] (if particle number is conserved)
            %3) spin sector [2*N_site + 1] (if spin is conserved (indepedendent of conserved particle number)
            %4) acting site of creation / annihilation operator
            %5) acting spin of -""- [2] (if not spinless)
            if obj.Basis.N_spin == 1
                
                operator_spin_range = 1;
                if ismember('N', obj.Structure)
                    particle_sector_range = obj.Basis.N_site + 1;
                else
                    particle_sector_range = 1;
                end
                
            else
                operator_spin_range = 2;
                if ismember('N', obj.Structure)
                    particle_sector_range = 2*obj.Basis.N_site + 1;
                else
                    particle_sector_range = 1;
                end
                
            end
            
            
            for dagger = -1 %NOTE: postive dagger version created by adjoint operation (')
                for c_site = 1:obj.Basis.N_site
                    %NOTE: start with at least one particle on the right side for the annhiliation operator
                    for N_part_ind = 2:particle_sector_range
                        N_part = N_part_ind - 1; 
                        if ismember('S', obj.Structure)
                            spin_range =  obj.Basis.f_spin_range(N_part);
                            spin_index_range = 1:numel(spin_range);
                        else
                            spin_index_range = 1;
                        end
                        %max(spin_index_range)
                        for Spin_ind = spin_index_range
                            for c_spin = 1:operator_spin_range
                                if obj.Basis.N_spin == 1
                                    %NOTE: spinless system
                                    if ismember('N', obj.Structure)
                                        
                                        if N_part - 1 >= 0
                                            Qmat_temp = f_c(obj, dagger, N_part, [], c_site, [], energy_cut);
                                            obj.Qmatrix{1, N_part_ind, Spin_ind, c_site, c_spin} = Qmat_temp;
                                            obj.Qmatrix{2, N_part_ind-1, Spin_ind, c_site, c_spin} = Qmat_temp';
                                        end
                                    else
                                            Qmat_temp = f_c(obj, dagger, [], [], c_site, [], energy_cut);
                                            obj.Qmatrix{1, N_part_ind, Spin_ind, c_site, c_spin} = Qmat_temp;
                                            obj.Qmatrix{2, N_part_ind, Spin_ind, c_site, c_spin} = Qmat_temp';
                                    end
                                
                                else
                                    if ismember('N', obj.Structure) && ismember('S', obj.Structure)
                                        
                                        %NOTE: compare definition taken from g_Qmat:
                                        %spin_sector_index = (Spin_inp - min(obj.Basis.f_spin_range(N_part_inp)))/2 + 1;
                                        
                                        Spin = 2*(Spin_ind - 1) + min(spin_range);
                                        new_Spin = Spin-2*(c_spin-1.5);
                                        new_Spin_ind = (new_Spin - min(obj.Basis.f_spin_range(N_part-1)))/2 + 1;
                                        
                                        if N_part - 1 >= 0 && ismember(new_Spin, obj.Basis.f_spin_range(N_part - 1))
                                            Qmat_temp = f_c(obj, dagger, N_part, Spin, c_site, 2*(c_spin-1.5), energy_cut);
                                            obj.Qmatrix{1, N_part_ind, Spin_ind, c_site, c_spin} = Qmat_temp;
                                            obj.Qmatrix{2, N_part_ind-1, new_Spin_ind, c_site, c_spin} = Qmat_temp';
                                        else
                                            warning(['new Spin: ', num2str(new_Spin), ' is out of range: ', num2str(obj.Basis.f_spin_range(N_part-1))])
                                        end
                                    elseif ismember('N', obj.Structure)
                                        if N_part - 1 >= 0
                                            Qmat_temp = f_c(obj, dagger, N_part, [], c_site, 2*(c_spin-1.5), energy_cut);
                                            
                                            obj.Qmatrix{1, N_part_ind, Spin_ind, c_site, c_spin} = Qmat_temp;
                                            obj.Qmatrix{2, N_part_ind-1, Spin_ind, c_site, c_spin} = Qmat_temp';
                                        end
                                    elseif ismember('S', obj.Structure)
                                        %check if Basis.f_spin_range() works without particle number
                                        Spin = 2*(Spin_ind - 1) + min(spin_range);
                                        new_Spin = Spin-2*(c_spin-1.5);
                                        new_Spin_ind = (new_Spin - min(spin_range))/2 + 1;
                                        if ismember(new_Spin, spin_range)
                                            Q_mat_temp = f_c(obj, dagger, [], Spin, c_site, 2*(c_spin-1.5), energy_cut);
                                            obj.Qmatrix{1, N_part_ind, Spin_ind, c_site, c_spin} = Q_mat_temp;
                                            obj.Qmatrix{1, N_part_ind, new_Spin_ind, c_site, c_spin} = Q_mat_temp';
                                        end
                                    else
                                        Qmat_temp = f_c(obj, dagger, [], [], c_site, 2*(c_spin-1.5), energy_cut);
                                            
                                        obj.Qmatrix{1, N_part_ind, Spin_ind, c_site, c_spin} = Qmat_temp;
                                        obj.Qmatrix{2, N_part_ind, Spin_ind, c_site, c_spin} = Qmat_temp';
                                    end
                                    
                                
                                end
                               
                                
                            end
                        end
                    end
                end
            end
            %{

            if ismember('S', obj.Structure), spin_conserved = true; else, spin_conserved = false; end
            if obj.Basis.N_spin == 1, spinless = true; else, spinless = false; end
            in_range = 0;
            not_in_range = 0;
            obj.Qmatrix = [];
            for dagger = -1 %NOTE: positive dagger version created by daggering
                for N_part_inp = 0:obj.Basis.N_spin*obj.Basis.N_site
                    
                    if spin_conserved && ~spinless
                        spin_range = obj.Basis.f_spin_range(N_part_inp);
                        spin_range2 = obj.Basis.f_spin_range(N_part_inp - 1); %dagger
                    else
                        spin_range = 0;
                        spin_range2 = 0; %dagger version
                        
                    end
                    if spinless
                        c_spin_range = -1;
                    else
                        c_spin_range = [-1,1];
                    end
                    
                    min_spin = min(spin_range);
                    min_spin2 = min(spin_range2); %dagger version
                    for Spin_index = (spin_range - min_spin + 1)
                        if obj.DBUG, disp(['Qmat generate, Dag: ', num2str(dagger), ', N: ', num2str(N_part_inp), ', S: ', num2str(Spin_index)]); end
                        for c_site = 1:obj.Basis.N_site
                            for c_spin_index = (c_spin_range/2 + 1.5)
                                
                                
                                
                                if spinless %spinless
                                    c_spin = [];
                                else
                                    c_spin = (c_spin_index-1.5) * 2;
                                end
                                
                                if spin_conserved && ~spinless %two spins, spin conserved
                                    Spin_ind = Spin_index + min_spin -1;
                                    Spin_index_dagger = Spin_ind - 1*c_spin - min_spin2 +1;
                                else
                                    Spin_ind = [];
                                    Spin_index_dagger = 1;
                                end
                                
                                
                                spin_range_dagger = obj.Basis.f_spin_range(N_part_inp + dagger);
                                if (N_part_inp + dagger) < 0 || (N_part_inp + dagger) >= obj.Basis.N_spin*obj.Basis.N_site || ...
                                        (~isempty(Spin_ind) && ~ismember(Spin_ind + dagger*c_spin, spin_range_dagger) )
                                    disp('not in range')
                                    not_in_range = not_in_range + 1;
                                else
                                    obj.Qmatrix{dagger/2+1.5, N_part_inp + 1, Spin_index, c_site, c_spin_index} = ...
                                        f_c(obj, dagger, N_part_inp, Spin_ind, c_site, c_spin, energy_cut);
                                     %dagger = +1
                                    if isempty(obj.Qmatrix{dagger/2+1.5, N_part_inp + 1, Spin_index, c_site, c_spin_index} )
                                        warning('empty Qmatrix entry!!!')
                                    end
                                    obj.Qmatrix{2, N_part_inp - 1 + 1, Spin_index_dagger, c_site, c_spin_index} = obj.Qmatrix{dagger/2+1.5, N_part_inp + 1, Spin_index, c_site, c_spin_index}' ;
              
                                    
                                    in_range = in_range + 1;
                                end
                            end
                        end
                    end
                end
            end
            %}
            
            %disp(['in range: ', num2str(in_range)])
            %disp(['not in range: ', num2str(not_in_range)])
%             % ----------- separated version
%             obj.Qmatrix = [];
%             % for comparison with old program, ignores Spin and creates Q-matrices for whole spin sectors
%             % makes sense if spin is not conserved
%             
%             for dagger = [-1,+1]
%                 for N_part_inp = 0:obj.Basis.N_spin*obj.Basis.N_site
%                     if spin_less % spin less
%                         
%                         for c_site = 1:obj.Basis.N_site
%                             obj.Qmatrix{dagger/2+1.5, N_part_inp + 1, 1, c_site, 1} = ...
%                                 f_c(obj, dagger, N_part_inp, [], c_site,[]);
%                         end
%                         
%                     elseif ~spin_conserved  % not spin conserved
%                         
%                         
%                         for c_site = 1:obj.Basis.N_site
%                             for c_spin = [-1,1]
%                                 obj.Qmatrix{dagger/2+1.5, N_part_inp + 1, 1, c_site, c_spin/2 + 1.5} = ...
%                                     f_c(obj, dagger, N_part_inp, [], c_site, c_spin);
%                             end
%                         end
%                     else
%                         
%                         spin_range = obj.Basis.f_spin_range(N_part_inp);
%                         min_spin = min(spin_range);
%                         for Spin_inp = spin_range
%                             for c_site = 1:obj.Basis.N_site
%                                 for c_spin = [-1,1]
%                                     obj.Qmatrix{dagger/2+1.5, N_part_inp + 1, Spin_inp-min_spin + 1, c_site, c_spin/2 + 1.5} = ...
%                                         f_c(obj, dagger, N_part_inp, Spin_inp, c_site, c_spin);
%                                 end 
%                             end  
%                         end
%                     end
%                 end 
%             end
        end
        
        function obj = f_reset_Qmatrix(obj)
            obj.Qmatrix = cell(2, 2*obj.Basis.N_site + 1, 2*obj.Basis.N_site+2, obj.Basis.N_site, 2);
            obj.Qmatrix_counter = 0;
        end
        
        %works for Spin and Spinless, quite new - used for local storage in cell (fast), also used in SCB
        function [Q_matrix, obj] = g_Qmat(obj, dagger, N_part_inp, Spin_inp, c_site, c_spin, energy_cut)
            if nargin < 7, energy_cut = []; end
            if ~iscell(obj.Qmatrix)
                error('perform Energy.f_Qmatrix_generate first') % or initialize in constructor of Energy
            end
            %spinless: Spin_inp & c_spin empty
            %spin not conserved Spin_inp empty
            %particle number not conserved N_part_inp empty
            
            %NOTE: some warnings to check if use of g_Qmat is done correctly 
            %(check if particle / spin sector arguments match conserved quantities):
            if obj.Basis.N_spin == 1
                if ~isempty(Spin_inp) || ~isempty(c_spin)
                    warning('In spinless case g_Qmat should not be called with Spin arguments!')
                end
            else
                if ~ismember('N', obj.Structure) && ~isempty(N_part_inp)
                    warning('If particle number is not conserved Q-matrices for a special particle sector are not meaningful.')
                elseif ~ismember('S', obj.Structure) && ~isempty(Spin_inp)
                    warning('If spin number is not conserved Q-matrices for a special spin sector are not meaningful.')
                end
            end
            
            if ismember('N', obj.Structure) && isempty(N_part_inp)
                error('preserved particle number has to be set, use g_full_Qmat for return of full Qmatrix!');
            end
            if ismember('S', obj.Structure) && isempty(Spin_inp)
                error('preserved spin number has to be set, use g_full_Qmat for return of full Qmatrix!');
            end
            
            
            %NOTE: Qmatrix indices:
            %1) no dagger, dagger [2]
            %2) particle_sector_index [basis.N_spin*basis.N_site+1] (if particle number is conserved)
            %3) spin_sector_index [2*N_site + 1] (if spin is conserved (indepedendent of conserved particle number)
            %4) operator_site_index of creation / annihilation operator
            %5) operator_spin_index -""- [2] (if not spinless)
            
            %NOTE: possibilities for Qmatrix indices:
            
            %particle_sector_index = N_part_inp + 1;
            %particle_sector_index = 1;
            
            %spin_sector_index = Spin_inp - min(obj.Basis.f_spin_range(N_part_inp)) + 1;
            %spin_sector_index = 1;
            
            %operator_spin_index = c_spin/2+1.5;
            %operator_spin_index = 1;
            
            if obj.Basis.N_spin == 1
                spin_sector_index = 1;
                operator_spin_index = 1;
                if ismember('N', obj.Structure)
                    particle_sector_index = N_part_inp + 1; 
                else
                    particle_sector_index = 1;
                end
            else
                operator_spin_index = c_spin/2+1.5;
                if ismember('N', obj.Structure)
                    particle_sector_index = N_part_inp + 1; 
                else
                    particle_sector_index = 1;
                end
                if ismember('S', obj.Structure)
                    spin_sector_index = (Spin_inp - min(obj.Basis.f_spin_range(N_part_inp)))/2 + 1;
                else
                    spin_sector_index = 1;
                end
            end
            
            
            %NOTE: works for all cases. f_c used as internal method, no checking in that function, since check already performed in
            % g_Qmat
            %TODO: make f_c a real private method.
            
            
            
            Q_matrix = obj.Qmatrix{dagger/2+1.5, particle_sector_index, spin_sector_index, c_site, operator_spin_index};
            if isempty(Q_matrix)
                
                Q_matrix = obj.f_c(dagger, N_part_inp, Spin_inp, c_site, c_spin, energy_cut);
                obj.Qmatrix{dagger/2+1.5, particle_sector_index, spin_sector_index, c_site, operator_spin_index}  =  Q_matrix;
            end

        end
        
        
        % returns full Q-matrix
        function Q_matrix = g_full_Qmat(obj,dagger, c_site,c_spin)
            Q_matrix = sparse (obj.Numel_En, obj.Numel_En);
            %NOTE: implement different setups (spinless, spin/particle conserved)
            if obj.Basis.N_spin > 1
                if ismember('N', obj.Structure) && ismember('S', obj.Structure)
                    for n_part = (0:obj.Basis.N_site*obj.Basis.N_spin-1)  + (1-dagger)/2
                        for spin = obj.Basis.f_spin_range(n_part)
                            
                            right_index = obj.f_index('NS',n_part, spin);
                            
                            left_index = obj.f_index('NS',n_part + dagger, spin + dagger*c_spin);
                            if ~isempty(left_index)
                                Q_matrix(left_index,right_index) = obj.g_Qmat(dagger, n_part, spin, c_site, c_spin);
                            end
                        end
                    end
                elseif ismember('N', obj.Structure)
                    for n_part = (0:obj.Basis.N_site*obj.Basis.N_spin-1)  + (1-dagger)/2
                        right_index = obj.f_index('N',n_part);
                        left_index = obj.f_index('N',n_part + dagger);
                        if ~isempty(left_index)
                            Q_matrix(left_index,right_index) = obj.g_Qmat(dagger, n_part, [], c_site, c_spin);
                        end
                    end
                        
                elseif ismember('S', obj.Structure)
                    for spin = obj.Basis.f_spin_range()
                        right_index = obj.f_index('S',spin);
                        left_index = obj.f_index('S',spin + dagger*c_spin);
                        if ~isemtpy(left_index)
                            Q_matrix(left_index,right_index) = obj.g_Qmat(dagger, [], spin, c_site, c_spin);
                        end
                    end
                else
                    right_index = obj.f_index();
                    left_index = obj.f_index();
                    if ~isemtpy(left_index)
                        Q_matrix(left_index,right_index) = obj.g_Qmat(dagger, [], [], c_site, c_spin);
                    end
                    
                end
            else
                if ismember('N', obj.Structure)
                    for n_part = (0:obj.Basis.N_site-1)  + (1-dagger)/2
                        right_index = obj.f_index('N',n_part);
                        left_index = obj.f_index('N',n_part + dagger);
                        if ~isempty(left_index)
                            Q_matrix(left_index,right_index) = obj.g_Qmat(dagger, n_part, [], c_site, []);
                        end
                    end
                    
                    
                 else
                    right_index = obj.f_index();
                    left_index = obj.f_index();
                    if ~isempty(left_index)
                        Q_matrix(left_index,right_index) = obj.g_Qmat(dagger, [], [], c_site, []);
                    end
                    
                end
            end
        end
        
        
        function Q_all = g_all_full_Qmat(obj)
            %TODO: rewrite!
            % make a similar structure as for Q_matrix storage!
            % here full Q_matrices are stored, 
            % maybe this function is used for self consistent Born?
            Q_all = cell(obj.Basis.N_site, 2, obj.Basis.N_spin);
            for dagger = -1:2:2
                for c_site = 1:obj.Basis.N_site
                    for c_spin = -1:2:(obj.Basis.N_spin -1)
                        Q_all{c_site,dagger/2+1.5,c_spin/2+1.5} = obj.g_full_Qmat(dagger, c_site, c_spin);
                    end
                    
                end
                
            end
        end
        
        function [sub_index] = f_energy_subindex(obj, n_part, spin, energ_range)
            %TODO: handle empty spin
            index_all = obj.Index('NS', n_part, spin);
            if ismember('A', obj.Structure)
                tag = 'NSA'; 
            else
                tag = 'NSE';
            end
            
            index_energies = obj.Index(tag, n_part, spin, energ_range);
            warning('Replace by Subindex function')
            %TODO: if this takes tooo long replace by a different scheme:
            % subtract minimum index number of block - should have the same result
            [~,sub_index] = ismember(index_energies, index_all);
     
        end
        
        function [sub_index, sub_index_2] = Subindex(obj, First_filter, First_arguments, Second_filter, Second_arguments)
            %Subindex: returns the subindices identified by the second filter in all elements identified by the
            % e.g. Subindex('NS', {npart, spin}, 'E', {{-inf, E_cut}})
            %first filter
            
            small_index = obj.Index([First_filter, Second_filter], First_arguments{:}, Second_arguments{:});
            large_index = obj.Index(First_filter, First_arguments{:});
            [~,sub_index] = ismember(small_index, large_index);
            
            %Note: alternative!
            % if Categ is ordered in the same order use the following function instead:
            % (avoids the expensive ismember function)
            sub_index_2 = small_index - min(large_index)+1;
            
            
        end
        

        % conversion to system sites
        %TODO: old - rewrite!!!!
        function[tab] = s_site_info(obj, varargin )
            if nargin < 2 || isempty(varargin)
                index = true(obj.Numel_En,1);
            else
                
                index = obj.Categ.filter_inp(varargin{:}); %find(filter_inp(varargin(:), obj.Categ));
            end
            Table_label = {'Energ', obj.Basis.Binary_Label{:}};
            if ismember('S', obj.Structure)
                Spin_column = obj.Spin(index);
                Table_label = {'Spin', Table_label{:}};
            else
                Spin_column = zeros(numel(index),0);
            end
            if ismember('N', obj.Structure)
                N_part_column = obj.N_part(index);
                Table_label = {'N_part', Table_label{:}};
            else
                N_part_column = zeros(numel(index),0);
            end
            
            
            tab = array2table([N_part_column, Spin_column , obj.Energies(index), ((obj.Eigenvector(:,index).^2).')*double(obj.Basis.f_bin())], 'VariableNames', Table_label );
            
            
%             % site_info  - method - gives representation of certain Eigenvectors in spatial space representation
%             %   site_info(obj, En_input, Np_input, Sp_input )
%             %
%             if nargin < 3 || isempty(Np_input), L_N = true;
%             else, L_N = ismember(obj.N_part, Np_input); end  %if ismember is too slow use internal function ismembc (in combination with sort)
%             if nargin < 4 || isempty(Sp_input), L_S = true;
%             else, L_S = ismember(obj.Spin,Sp_input); end
%             if nargin < 2 || isempty(En_input), L_E = true;
%             else, L_E = ismember(obj.Energies,En_input); end
%             L = L_N & L_S & L_E & true(size(obj.Energies));
%             
%             tab = array2table(((obj.Eigenvector(:,L).^2).'*obj.Basis.f_bin()),'VariableNames', obj.Basis.Binary_Label);
%             %disp(array2table(((obj.Eigenvector(:,L).^2).'*obj.Basis.f_bin),'VariableNames', obj.Basis.Binary_Label))
        end
        
        
        
        %TODO: old - rewrite
        function obj = f_cut_Ener(obj,varargin)
            index = find(filter_inp(varargin(:), obj.Categ));
            warning('check if outdated - some properties went to Categories!!!')
            
            obj.Energies = obj.Energies(index);
            obj.N_part = obj.N_part(index); 
            obj.Spin = obj.Spin(index);
            obj.Eigenvector = obj.Eigenvector(:,index);
            obj.Index = obj.Index(index,:);
            obj.N_part_var = obj.N_part_var(index);
            obj.Spin_var = obj.Spin_var(index);
            obj.Degeneracy = obj.Degeneracy(index);
            
            for n_Cat = 1: numel(obj.Categ)
                obj.Categ(n_Cat).Value = obj.Categ(n_Cat).Value(index);
            end
            
            
        end
        
        

        
        % return special energies, depending on N,S, etc...
        function Energ = g_En(obj, Np_input, Sp_input)
            % g_En(obj, Np_input, Sp_input)
            %
            % En  - getter - returns Energy
            % Np_input and  Sp_input work as scalar filters. if not set g_En returns all Energies
            % 
            % if Energy has no Spin property (spinless case), spin is ignored.
            %NOTE: new and checked :)
            warning('Outdated')
            if nargin < 2 || isempty(Np_input),  L_N = true; else,  L_N = obj.N_part == Np_input; end
            if nargin < 3 || isempty(Sp_input)
                L_S = true; 
            else
                if obj.basis.N_spin > 1
                    L_S = obj.Spin == Sp_input; 
                else
                    L_S = true;
                end
            end

            L = L_N & L_S & true(size(obj.Energies));
            
            Energ = obj.Energies(L);
        end
        
        
        function [Eigenvec] = f_Ev(obj, Np_input, Sp_input, Deg_input)
            % [Eigenvec] = f_Ev(obj, Np_input, Sp_input, Deg_input)
            %
            % Ev  - method - returns the eigenvectors according to a special part.sec./Spin/Degeneracy
            % 
            %
            % if Energy has no Spin property (spinless case), spin is ignored.
            
            error('maybe old');
            
            
            if nargin < 2 || isempty(Np_input), L_N = true;
            else, L_N = ismember(obj.N_part, Np_input); end  %if ismember is too slow use internal function ismembc (in combination with sort)
            if nargin < 3 || isempty(Sp_input), L_S = true;
            else, L_S = ismember(obj.Spin,Sp_input); end
            if nargin < 4 || isempty(Deg_input), L_D = true;
            else, L_D = ismember(obj.Degeneracy,Deg_input); end
            L = L_N & L_S & L_D & true(size(obj.Energies));
            
            L2 = obj.Basis.f_ind(Np_input, Sp_input);
            
            Eigenvec = obj.Eigenvector(L2,L);
            
        end
        
        function [Eigenvec] = g_Ev(obj, Np_input, Sp_input) 
            % [Eigenvec] = g_Ev(obj, Np_input, Sp_input)
            %
            % Ev  - method - returns the eigenvectors according to a special part.sec./Spin/Degeneracy
            %   
            % if Energy has no Spin property (spinless case), spin is ignored.
            %
            % TODO: add by sector, write a summary of possibilities to address Eigenvectors (either full eigenvector or just rows living in sector)
            % NOTE: new and checked :)
            % compare with Eigen_Vector - Think about it!!!
              warning('Outdated')
            
            if nargin < 2 || isempty(Np_input),  L_N = true; else,  L_N = obj.N_part == Np_input; end
            if nargin < 3 || isempty(Sp_input)
                L_S = true; 
            else
                if obj.Basis.N_spin > 1
                    L_S = obj.Spin == Sp_input; 
                else
                    L_S = true;
                end
            end
            
            L = L_N & L_S & true(size(obj.Energies));
            
            %TODO: Big question mark - how does this work - how to put this in one framework, deal with spinless
            %systems - think about it!!!
            %NOTE: now just changed the following line below. maybe a more logic approach can be done ;)
           
            L2 = obj.Basis.Index('NS', Np_input, Sp_input);
            
            Eigenvec = obj.Eigenvector(L2,L);
        end
        
        
        function [index] = f_In(obj, Np_input, Sp_input, Deg_input)
            % In  - method - returns the Index of the Eigenstates according to a special part.sec./Spin/Degeneracy
            %   In(obj, Np_input, Sp_input, Deg_input)
            %
            if nargin < 2 || isempty(Np_input), L_N = true;
            else, L_N = ismember(obj.N_part, Np_input); end  %if ismember is too slow use internal function ismembc (in combination with sort)
            if nargin < 3 || isempty(Sp_input), L_S = true;
            else, L_S = ismember(obj.Spin,Sp_input); end
            if nargin < 4 || isempty(Deg_input), L_D = true;
            else, L_D = ismember(obj.Degeneracy,Deg_input); end
            L = L_N & L_S & L_D & true(size(obj.Energies));
            
            index = obj.Index(L);
            
            error('maybe old')
        end
        
        function [index] = f_index(obj,varargin)
            %NOTE: changed from external filter_inp function to Categ.filter_inp function
            if nargin == 1
                index = obj.Index;
            else
                L = obj.Categ.filter_inp(varargin{:}); 
                index = find(L);
            end
            
        end
        
        function [Tab] = s_print(obj, varargin)
            % print  - method - prints a table of the current Energies with part. sec / Spin / Deg attributes
            %    print(option)
            %    the option can be chosen to be 'l' to display the lowest energies and 'a' to print all energies
            if nargin < 2  || isempty(varargin), varargin = {'a'}; end
            
            if numel(varargin) == 1
                option = varargin{1};
            
                switch lower(option(1))
                    case 'a' % show all energies
                        if ismember('S', obj.Structure) && ismember('N', obj.Structure)
                            Tab = table(obj.Energies, obj.N_part, obj.Spin, obj.Degeneracy, 'variablenames',{'Energies', 'Part_numb', 'Spin', 'Degeneracy'});
                        elseif ismember('S', obj.Structure)
                            Tab = table(obj.Energies, obj.Spin, obj.Degeneracy, 'variablenames',{'Energies', 'Spin', 'Degeneracy'});
                        elseif ismember('N', obj.Structure)
                            Tab = table(obj.Energies, obj.N_part, obj.Degeneracy, 'variablenames',{'Energies', 'Part_numb', 'Degeneracy'});
                        else
                            Tab = table(obj.Energies, obj.Degeneracy, 'variablenames',{'Energies', 'Degeneracy'});
                        end
                    case 'l' % lowest 10 energies
                        % question whether to show unique 10 lowest energies
                        % and there degeneracy according to particle / spin sector
                        [~, ix]  = sort(obj.Energies);
                        index = 1:min(10,numel(obj.Energies));
                        if ismember('S', obj.Structure) && ismember('N', obj.Structure)
                            Tab = table(obj.Energies(ix(index)), obj.N_part(ix(index)), obj.Spin(ix(index)), obj.Degeneracy(ix(index)), 'variablenames',{'Energies', 'Part_numb', 'Spin', 'Degeneracy'});
                        elseif ismember('S', obj.Structure)
                            Tab = table(obj.Energies(ix(index)), obj.Spin(ix(index)), obj.Degeneracy(ix(index)), 'variablenames',{'Energies', 'Spin', 'Degeneracy'});
                        elseif ismember('N', obj.Structure)
                            Tab = table(obj.Energies(ix(index)), obj.N_part(ix(index)), obj.Degeneracy(ix(index)), 'variablenames',{'Energies', 'Part_numb', 'Degeneracy'});
                        else
                            Tab = table(obj.Energies(ix(index)), obj.Degeneracy(ix(index)), 'variablenames',{'Energies', 'Degeneracy'});
                        end
                end
            else
                index = obj.Index(varargin{:});
                
                Tab_array = zeros(numel(index), obj.Categ.Numel_Data_Types);
                Tab_name = cell(1, obj.Categ.Numel_Data_Types);
                for k = 1:obj.Categ.Numel_Data_Types
                    Tab_array(:,k) = obj.Categ.Data(k).Values(index);
                    Tab_name{1,k} = obj.Categ.Data(k).Name;
                end
                Tab = array2table(Tab_array, 'variablenames', Tab_name);
                %Tab = table(obj.Energies(index), obj.N_part(index), obj.Spin(index), obj.Degeneracy(index), 'variablenames',{'Energies', 'Part_numb', 'Spin', 'Degeneracy'});
            end
                
        end
        
        function [] = s_spectrum(obj,num_bins, total_degeneracy, details, ax)
            if nargin < 2 || isempty(num_bins), num_bins = 150; end
            if nargin < 3 || isempty(total_degeneracy), total_degeneracy = false; end
            if nargin < 4 || isempty(details), details = false; end
            if nargin < 5 || isempty(ax), new_figure = true; else new_figure = false; end
            
            if details == true || obj.Basis.N_site < 5
                if new_figure
                    figure()
                    plot_settings
                    ax = axes;
                else
                    axes(ax)
                end
                [sort_ener, ix] = sort(obj.Energies);
                if numel(sort_ener) > 300
                    index_lim = obj.Index('E', {min(obj.Energies)-10^-10, sort_ener(300)});
                else
                    index_lim = 1:numel(sort_ener);
                end
                
                n_part = obj.N_part(index_lim);
                energ = obj.Energies(index_lim);
                %TODO: rewrite for spinless systems!!!
                if ismember('S', obj.Structure)
                    spin = obj.Spin(index_lim);
                    max_spin = max(abs(spin));
                    spin_text = {'spin 0', 'spin 1', 'spin 2', 'spin 3', 'spin 4', 'spin 5', 'spin 6'};
                else
                    spin = [];
                    max_spin = 1;
                    spin_text = {};
                end
                
                tol_round = 10^-5;
                
                if total_degeneracy
                    degen = obj.Degeneracy;
                else
                    
                    X = [n_part, spin, round(energ/tol_round)*tol_round];
                    %X(2:3,:) = X([1,1],:)
                    %p = randperm(64)
                    %[A, B, C]= unique(X(p,:), 'rows')
                    [~,B, C] = unique(X, 'rows');
                    D = unique(C);
                    E = histc(C,D);
                    
                    if diff(B)< 0
                        warning('Energy not well sorted!')
                    end
                    %TODO: shrinks degen, doesn't duplicate the numbers for the same states - procede as done for
                    % obj.Degeneracy or apply for loops.
                    degen = E(C);
                    %                    [ener_unique,B,ind_position] = unique(X, 'rows');
                    %                     frequency = histc(X(:,3),ener_unique(:,3));
                    %                     degen = frequency(ind_position);
                    
                end
                
                
                degen(degen > 12) = 13;
                marker = {'x','<','^','s','p','h', '>', 'o', 'v', 'd', '+', '.', '*'};
                color_spin = colormap(gca, lines(7));
                
                xlimit = ax.XLim;
                xmode = ax.XLimMode;
                ylimit = ax.YLim;
                ymode = ax.YLimMode;
                index = 1:numel(index_lim);
                index(isinf(obj.Energies(index_lim) )) = [];
                for k_in = 1:numel(index)
                    k = index(k_in);
                    if ismember('S', obj.Structure)
                        spin_mark = spin(k);
                    else
                        spin_mark = 0;
                    end
                    
                    plot(n_part(k) + 0.3*spin_mark/max_spin, energ(k), marker{degen(k)}, 'color', color_spin(abs(spin_mark)+1,:))
                    if strcmp( xmode, 'manual')
                        xlim(xlimit)
                    end
                    if strcmp(ymode, 'manual')
                        ylim(ylimit)
                    end
                    hold on
                end
                
                
                for k=1:13
                    pl(k) = plot(nan,nan, marker{k}, 'color', 'k');
                    hold on
                end
                if ismember('S', obj.Structure)
                    for k=1:7
                        pl(13+k) = plot(nan,nan,'-', 'color', color_spin(k,:));
                    end
                end
                
                hold off
                legend(pl, '1x','2x','3x', '4x', '5x', '6x', '7x', '8x', '9x', '10x', '11x', '12x', 'more', spin_text{:}, 'location','northeastoutside')
                ylabel('Energy')
                xlabel('Particle Sector')
               %legend(pl(14:end), '$|n_\uparrow - n_\downarrow| = 0$', '$|n_\uparrow - n_\downarrow| = 1$', '$|n_\uparrow - n_\downarrow| = 2$', '$|n_\uparrow - n_\downarrow| = 3$', '$|n_\uparrow - n_\downarrow| = 4$', '$|n_\uparrow - n_\downarrow| = 5$', '$|n_\uparrow - n_\downarrow| = 6$', 'interpreter', 'latex', 'fontsize', 12)
                set(gcf,'color','w')
                grid on
                title('Energy spectrum of system Hamiltonian (Energy class)')
            else
                if new_figure
                    figure('position', [300,300,700,300])
                else
                    axes(ax)
                end
                energ_vec =  linspace(min(obj.Energies), max(obj.Energies(~isinf(obj.Energies))),num_bins);
                Bar = zeros(num_bins-1, obj.Basis.N_spin*obj.Basis.N_site+1);
                for k = 0:obj.Basis.N_spin*obj.Basis.N_site
                    energ = obj.Energies('N', k);
                    if numel(energ) > 2
                        Bar(:,k) = histcounts(energ,energ_vec );
                    end
                end
                bb = bar3(energ_vec(1:end-1)', Bar);
                
                for k = 1:length(bb)
                    zdata = bb(k).ZData;
                    bb(k).CData = zdata;
                    bb(k).FaceColor = 'interp';
                end
                cmap = colormap(parula(100));
                cmap(1,:) = [1,1,1];
                colormap(cmap)
                caxis([0, max(max(Bar(1:end-1,:)))])
                colorbar('southoutside')
                
                ylabel('Energy')
                xlabel('Particle Sector')
                view(-90,90)
                set(gcf,'color','w')
                
            end
        end
    end
    methods(Static)
        function En = sf_from_system(en_input)
            basis = en_input.basis;
            system = en_input.system;
            basis = orderfields(basis, {'N_site', 'N_sector', 'S_sector', 'Docc_input', 'N_spin'});
            bas_param = struct2cell(basis);
            Bas = Basis(bas_param{:});
            
            if ~isfield(system, 'class')
                system.class = 'Hubbard';
                Geom = Geometry(Bas.N_site, 'full', ones(Bas.N_site));
                system.geometry = Geom;
                system.bas = Bas;
                system = orderfields(system, {'class', 'geometry', 'b', 'U', 'xi', 'V', 'bas', 'part_hole_sym'});%; fn(1:end-1)]);
            end
            %system = orderfields(system, [{'class', 'b', 'U', 'xi', 'V', 'part_hole_sym'}]);%; fn(1:end-1)]);
            sys_param = struct2cell(system);
            Sys = eval([sys_param{1}, '(sys_param{2:end})']);
            
            En = Sys.f_solve;
        end
    end
        
        % old functions
        
%         function [block_table,L] = f_block_table(obj, N_part_inp, Spin_inp, Ener_inp, order)
%             % [block_table,L] = f_block_table(obj, N_part_inp, Spin_inp, Ener_inp, order)
%             % 
%             % TODO: reevaluate if needed or not coped by Tab and Filter functions
%             % block_table - method - creates table of groups according to N, S, E
%             %   block_table(obj, N_part_inp, Spin_inp, Ener_inp, order)
%             %   N_part_inp, Spin_inp and Ener_inp used to filter, can be vectors, an empty vector [] considers all
%             %   order can be set to group eigenstates in sectors, also defines the order
%             %   block_table illustrates the selected sectors
%             %   max_block_size illustrates how many energies are in each group
%             %   fil_block_size illustrates how many energies are in each group according to given filters
%             
%             
%             
%             % TODO: introduce a range feature, maybe with a cell: {from , till}
%             
%             % TODO: check if this also works for a differently ordered energies ('NES', 'NE', 'EN')
%             
%             if nargin < 5 || isempty(order), order = 'NSE'; end
%             
%             if nargin < 4 || isempty(Ener_inp), L_E = true; Ener_inp = [];
%             else L_E = ismembertol(obj.En, Ener_inp, 10^-12); end
%             if nargin < 3 || isempty(Spin_inp), L_S = true; Spin_inp = [];
%             else L_S = ismember(obj.Spin, Spin_inp); end
%             if nargin < 2 || isempty(N_part_inp), L_N = true; N_part_inp = [];
%             else L_N = ismember(obj.N_part, N_part_inp); end
%             
%             [use_N, N_index] = ismember('N',order);
%             [use_E, E_index] = ismember('E',order);
%             [use_S, S_index] = ismember('S',order);
%             
%             
%             %sort order:
%             num_of_col = use_N+use_S+use_E;
%             if use_N, col_num(N_index) = 2; tab_col_name{N_index} = 'N_part'; end
%             if use_E, col_num(E_index) = 1; tab_col_name{E_index} = 'Energ'; end
%             if use_S, col_num(S_index) = 3; tab_col_name{S_index} = 'Spin'; end
%             tab_col_name{num_of_col+1} = 'max_block_size';
%             tab_col_name{num_of_col+2} = 'fil_block_size';
%             
%             if num_of_col ~= numel(col_num),  error('wrong order index, allowed is N,E or S');end
%             
%             
%             %TODO: make more generell (what happens when ordering of Energy class is different?)
%             Ener_tab = table2array(obj.s_print) ;
%             
%             [Block_table,~, frequ] = unique(Ener_tab(:,col_num),'rows');
%             
%             Block_table = [Block_table, histc(frequ,1:size(Block_table,1))];
%             filter_bl_col = size(Block_table,2) + 1; %column for storing actual block size (incl. filtering)
%             
%             % use Block_table to filter special blocks of N_part, Spin or Energy
%             % TODO: this function may be slow - implement an if-clause for not filtering situation
%             [LB_N, LB_E, LB_S] = deal(true);
%             if use_N && ~isempty(N_part_inp) && all(isnumeric(N_part_inp)) , LB_N = ismember(Block_table(:,N_index), N_part_inp); end
%             if use_S && ~isempty(Spin_inp) && all(isnumeric(Spin_inp)), LB_S = ismember(Block_table(:,S_index), Spin_inp); end
%             if use_E && ~isempty(Ener_inp) && all(isnumeric(Ener_inp)), LB_E = ismember(Block_table(:,E_index), Ener_inp); end
%             LB = LB_N & LB_E & LB_S & true(size(Block_table,1),1);
%             Block_table = Block_table(LB,:);  %restrict to wished blocks
%             
%             % TODO: implement tolerance?
%             % TODO: implement here the indexing scheme and save those in an extra column
%             % TODO: implement diagonal elements
%             % TODO: think of how to transfere back (structure and filter)
%             L = cell(sum(LB),1);
%             index = 0;
%             for n_block = 1:sum(LB)
%                 if use_N, L_N = Block_table(n_block,N_index) == obj.N_part; end
%                 if use_E, L_E = abs(Block_table(n_block,E_index)-obj.Energies) < 10^-12; end
%                 if use_S, L_S = Block_table(n_block,S_index) == obj.Spin;  end
%                 L{n_block} = L_N & L_E & L_S;
%                 block_size = sum(L{n_block});
%                 Block_table(n_block, filter_bl_col) = block_size;
%                 Block_table_cell{n_block,1} = index + (1:block_size^2); %all block indizes (normal matlab indexing)
%                 Block_table_cell{n_block,2} = index + (1:(block_size+1):block_size^2); %diagonal terms (either as index or as logical)
%                 index = index + sum(L{n_block})^2;
%                 
%             end
%             
%             
%             block_table = array2table(Block_table, 'variablenames', tab_col_name);
%             [Block_table_cell{:,2}] %diagonal terms
%             
%             
%         end
        
        

        %         % return special energies, depending on N,S, etc...
%         % TODO: maybe old - see .Energies (via Categories
%         function [Energ, Part, Spi, Deg, Eigenvec] = f_En(obj, Np_input, Sp_input, Deg_input)
%             % [Energ, Part, Spi, Deg, Eigenvec] = f_En(obj, Np_input, Sp_input, Deg_input)
%             % 
%             % En  - method - returns Energy, particle sector, Spin, Degeneracy, and the eigenvectors
%             % TODO: apply filter method
%             %
%             if nargin < 2 || isempty(Np_input), L_N = true;
%             else
%                 if numel(Np_input) == 1
%                     L_N = obj.N_part == Np_input;
%                 else
%                     L_N = ismember(obj.N_part, Np_input);
%                 end  %if ismember is too slow use internal function ismembc (in combination with sort)
%             end
%             if nargin < 3 || isempty(Sp_input), L_S = true;
%             else
%                 if numel(Sp_input) == 1
%                     L_S = obj.Spin == Sp_input;
%                 else
%                     L_S = ismember(obj.Spin,Sp_input);
%                 end
%             end
%             if nargin < 4 || isempty(Deg_input), L_D = true;
%             else
%                 if numel(Deg_input) == 1
%                     L_D = obj.Degeneracy == Deg_input;
%                 else
%                     L_D = ismember(obj.Degeneracy,Deg_input);
%                 end
%             end
%             L = L_N & L_S & L_D & true(size(obj.Energies));
%             
%             Energ = obj.Energies(L);
%             Part = obj.N_part(L);
%             Spi = obj.Spin(L);
%             Deg = obj.Degeneracy(L);
%             Eigenvec = obj.Eigenvector(:,L);
%             
%         end
        
                %NOTE: OLD: creates a Qtable (to lookup), there is a new class Q_Table < Table (necessary?) maybe used in SCB
%         function Q_table = f_create_Q_table(obj, energy_cut)
%             error('OLD FUNCTION');
%             if nargin < 2, energy_cut = []; end
%             % TODO: old use class Q_Table instead, does the same?
%             % f_create_Q_table  - method - creates Q_table
%             %   create_Q_table()
%             %
%             N_ket_final = [];
%             Npart = [];
%             spin = [];
%             N_bra_final = cell(1,4);
%             for Np = 0:obj.Basis.N_spin*obj.Basis.N_site
%                 Spi = obj.Basis.f_spin_range(Np)';
%                 
%                 N_ket = zeros(size(Spi));
%                 for dagger = [-1,1]
%                     for c_spin = [-1,1]
%                         q_ind = dagger+2 + c_spin/2+0.5;
%                         N_bra{q_ind} = zeros(size(Spi));
%                         
%                         for k = 1:numel(Spi)
%                             N_ket(k) = numel(obj.g_En(Np,Spi(k)));
%                             N_bra{q_ind}(k) = numel(obj.g_En(Np+dagger, Spi(k) + dagger*c_spin));
%                         end
%                         N_bra_final{q_ind} = [N_bra_final{q_ind}; N_bra{q_ind}];
%                     end
%                 end
%                 N_ket_final = [N_ket_final; N_ket];
%                 Npart = [Npart; Np*ones(size(Spi))];
%                 spin = [spin; Spi];
%             end
%             
%             %Q_table.N_ket_final = N_ket_final;
%             %Q_table.N_bra_final = N_bra_final;
%             
%             Q_table = Qtable(Npart, spin, N_ket_final, N_bra_final);
%         end
%         

%         function energ = f_Ener(obj,varargin)
%             error('USE obj.Energies')
%             %energ =  obj.Categ.f_access('Energies', varargin{:});
%         end

        

        %TO BE REPAIRED
%        function Q_index = Qindex(obj, dagger, c_spin, Np_input, Sp_input)
%            warning(' Qindex function - to be REPAIRED')
%            % Q_index  - method - returns a table to show how many non-zero Q-matrix elements each block has
%            %   Qindex(obj, dagger, c_spin, Np_input, Sp_input)
%             if nargin < 4 || isempty(Np_input) || nargin < 5 || isempty(Sp_input)
%                 Tab = [];
%                 for Np = 0:2*obj.Basis.N_site
%                     Spi = obj.Basis.f_spin_range(Np);
%                     
%                     N_ket = zeros(size(Spi));
%                     
%                     for k = 1:numel(Spi)
%                         N_ket(k) = numel(obj.En(Np,Spi(k)));
%                         num2(k) = numel(obj.En(Np + dagger, Spi(k) + dagger*c_spin));
%                     end
%                     
%                     
%                     Tab = [Tab; ones(numel(Spi),1)*Np, Spi(:),N_ket(:),num2(:), N_ket(:).*num2(:) ];
%                     
%                 end
%                 Q_index = array2table(Tab, 'variablenames', {'Np', 'Spin', 'numel_ket', 'numel_bra', 'total_elements'});
%                 
%             else
%                 
%                 
%             end
%        end
        

        
    
    
    
end