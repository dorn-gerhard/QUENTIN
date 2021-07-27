classdef Sigma
    
    %author: Gerhard Dorn
    %date: october 15
    %use: handle the reduced density operator sigma
    %features: 
    %          handle:
    %             - indexing scheme, that allows to operate in individual particle/spin/energy sector (according to
    %               Energy basis set)
    %             - indexing scheme, that handles superoperator acting on sigma - 
    %             - transforming of the two schemes
    %          analyse:
    %             - print / visualize Sigma - which sector is dominant - offer a position space representation
    %             - checks if block structure is available - if inter block states are populated significantly
    %               (address blocks of different particle nnmber / different Energy / different Spin)
    %             - degeneracy, coherence, entanglement
    %          initialize:
    %             - create a random valid density matrix with certains properties (block structure)
    %             - check sigma (hermitian, trace, positive definite) 
    %             - correct sigma (as steady state solution) (maybe when calculated)
    % general description of the system under influence of the baths, expression in terms of the energies
    
    properties
        Values              % full representation as sparse matrix in the Basis scheme
        Index               % maybe not necessary as lookup possible
        Tolerance = 10^-12; %
        Eigenvalues         % Eigenvalues of Sigma (when analysed)
        Eigenvectors        % Eigenvectors of Sigma (when analysed)
        Block_structure     % Represents the structure of Sigma (sorting / indexing)
        Analysis            % Values concerning Hermiticity, Trace and Positivity
        Tab                 % Table according to which Sigma is organized
        Energy_Cut          % defines if during Sigma calculation an energy cut was used
        Energy_Cut_2        % defines if during Sigma calculation a second energy cut (virtual excitations) was used
    end
    
    
    methods
        function obj = Sigma(Ener_inp, A, Tab_inp, Energy_Cut_inp, Energy_Cut_2_inp)
            % Sigma(Ener_inp, A, Tab_inp)
            %
            % DESC: 1) generate random matrix according to Tab, prioritize use of Tab, if none given take structure
            % of En - for testing issues
            %       2) structure unknown, handed over is A as matrix + Energy
            %       3) stucture known, A handed over as block or vector or matrix, Tab + Energy
            %
            % third argument optional, A can be already block or vector structure, then Tab must be set
            % otherwise A can be general square matrix which will then be analyzed
            % Energy not stored anymore (too big)
           
            % TODO: make connection to 
            % SAVE SIGMA AS BLOCKS INCLUDING TABLE
            % TODO: check whether storing as blocks or as sparse is faster / difference?
            
            if nargin < 4 || isempty(Energy_Cut_inp), Energy_Cut_inp = 10^10; end
            if nargin < 5 || isempty(Energy_Cut_2_inp), Energy_Cut_2_inp = Energy_Cut_inp; end
            
            
            % constructor
            if nargin < 2 || isempty(A)
                if nargin < 3 || isempty(Tab_inp)
                    Tab_inp = Table(Ener_inp, Ener_inp.Structure);
                end
                
                % -------------- randomly generated Sigma ----------------
                % create random valid density matrix
                % in Blockstructure
                % get Blocktable
                
                % defines structure of sigma, which sectors are used? what is the order? (first three options) what is the block
                % structure? (NSE)
                % TODO: enable different block sigma types for random generation? -> like Sigma(En,[],'NE')
                block_structure = Tab_inp.Order;
                obj.Block_structure = block_structure; %looks like to be the same as Energy.Structure!!!
                obj.Tab = Tab_inp;
                if ismember('E', block_structure)
                    energy_tag = 'E';
                else 
                    energy_tag = 'A';
                end
                %block_table = Ener_inp.f_block_table([],[],[], obj.Block_structure);
                
                % TODO: use block_table from f_analyse or from Table class
                % how to address info by block number?
                
                %A = obj.Tab.f_block(rand(obj.Tab.Numel_Vector_Elements ,1)-0.5 + (0.1i * rand(obj.Tab.Numel_Vector_Elements,1))); %real part dominant
                en_block_index = obj.Tab.Block_Index(energy_tag, {-4000000, Energy_Cut_inp});
                
                num_of_blocks = numel(en_block_index);
                block_weight = rand(num_of_blocks,1);
                block_weight = block_weight/sum(block_weight);
                
                for n_block_ind = 1:num_of_blocks
                    n_block = en_block_index(n_block_ind);
                    temp = rand(obj.Tab.Block_Size(n_block))-0.5 + 0.1i * rand(obj.Tab.Block_Size(n_block))-0.05;
                    
                    %rho = (A{n_block}' + A{n_block})/2; % hermitian
                    rho = (temp' + temp)/2;
                    [ev,ew] = eig(rho);
                    rho = ev*abs(ew)*ev';  %pos def
                    rho = rho/trace(rho) * block_weight(n_block_ind); %different weight for each block, trace = 1
                    
                    A{n_block,1} = round(rho,15);
                end
                
                obj.Values = A;
                
                analyse.positivity = true; 
                analyse.is_herm = true;
                analyse.trace_one = true;
                analyse.trace = 1;
                analyse.trace_positive = true;
                analyse.sum_neg_diag = 0;
                analyse.sum_neg_ev = 0;
                analyse.hermiticity_diff = 0;
                
                obj.Analysis = analyse;
                
                
                
            else
                % ------------- not randomly generated Sigma -------------------
                % procedure: check first if cell block or not
                % if it is cell block you need a Table, else error
                % if it is a matrix: is there a Table given
                % if not: analyse for block structure, use this one to store
                
                if iscell(A)
                    if nargin < 3 || isempty(Tab_inp), error('There has to be a defined Tablestructure for block_input'); end
                    block_structure = Tab_inp.Order;
                    obj.Values = A;
                    
                elseif ismatrix(A)
                    if any(size(A) == 1)
                        % A is a vector, check if Tab is set
                        if nargin < 3 || isempty(Tab_inp), error('There has to be a defined Tablestructure for vector_input'); end
                    
                    else
                        % A is really a matrix, check if Tab is set, otherwise analyse block structure
                        if nargin < 3 || isempty(Tab_inp)
                            
                            [block, block_structure, p2, block_column_label] = obj.f_analyse_block(A, Ener_inp);
                            disp(['Analysed block structure is: ', block_structure])
                            T_block =cell2table(block, 'variablenames', block_column_label) ;
                            disp(T_block)
                            %maybe save T_block as well?
                            
                            %TODO distinguish between: find out what structure is possible from the point of nnz in
                            %sigma, compared to: Is this structure compatible with given Energy Structure (especially
                            %the order)
                            % check whether block is already ordered (p2 is sorted)
                            
                            
                            
                            if all(ismember(block_structure, Ener_inp.Structure)) && all(p2(1:end-1) < p2(2:end))
                                Tab_inp = Table(Ener_inp, block_structure);
                            else
                                Tab_inp = Table(Ener_inp, '');
                            end
                            %define Tab
                            
                        end
                    end
                    
                    block_structure = Tab_inp.Order;
                    obj.Values = Tab_inp.f_block(A);
                    
                end
                obj.Energy_Cut = Energy_Cut_inp;
                obj.Energy_Cut_2 = Energy_Cut_2_inp;
                obj.Tab = Tab_inp;
                
                
               
                
                
                [obj,analyse] = obj.f_validate();
                
                
                
                % TODO:
                % restructure Energy???
                % a function like Energy.f_reorder(order, structure: EN, NES, E) to update the structure of the Energy
                % sometimes Sigma has a different structure (mixes energies), additional blocks due to symmetries,
                % preordering done already in System class for generating Energy class
                
                
                
                obj.Block_structure = block_structure;
                obj.Analysis = analyse;
                % make density matrix corrections if necessary, checks on hermiticity, positivity, trace
                
            end
            
            weight = cellfun(@(x) sum(abs(diag(full(x)))), obj.Values);
            xweight = cellfun(@(x) sum(sum(abs(full(x) - diag(diag(x))))), obj.Values);
            obj.Tab = obj.Tab.f_add_categ({'Weight', 'W'}, weight);
            obj.Tab = obj.Tab.f_add_categ({'XWeight', 'X'}, xweight);
            
            obj.Index = (1:numel(A)).';
            %             obj.Eigenvalues = [];
            %             obj.Eigenvectors = [];
            
          

        end
        
        %maybe difficult with name: returns sparse of whole sigma, or cuts out
        %subblocks (ignores block structure)
        function [S] = g_matrix(obj, varargin)
            % [S] = g_matrix(obj, varargin)
            % 
            % returns full Sigma matrix or filtered area (SUBBLOCK)
            % takes Category from Energy class! (since Sigma is arranged in that way)
            
            
            S = obj.Tab.f_matrix(obj.Values);
            
            if nargin > 1
                %Apply Filter
                % TODO: find a more efficient way? if input
                L = obj.Tab.List_Index(varargin{:});
                S = S(L, L);
                
              
            end
                
            
            
        end
        
        
        % TODO: could be renamed to g_blocks (search in block structure), 
        % 4 output options: Index, Blocks, sparse Matrix, Vector
        % use nargout to save time! definition of different input: filter or indices 
        % (indices already included in Tab (should be first argument) - 
        % in Tab the sequence should be adopted (not include empty fields in Category!)
        
        % for special version (superblock feature (search in subspace): write own function that return indices!
        
        %TODO: OLD - change and rethink!!!!!===================================??????????????????
        function [Indices, Blocks, Matrix, Vector] = g_superblocks(obj, varargin)
            % S = g_superblocks(obj, varargin)
            % 
            % returns the blocks that fulfil filter criteria in table as indices, blocks, sparse matrix or vector
            % searches in Table, additionally properties not addressed
            % for example get me all blocks, that inherit a special subproperty (in energy class)
            % TODO: for that we need a special function that searces also in the subproperties if filter options
            % are not addressed in Tab.Category (return-variable of function filter_inp containing rest filter)
            % 
            
            if nargin > 1
                
                L = obj.Tab.Block_Index(varargin{:});
                Indices = L;
                %L = filter_inp(varargin, obj.Tab.Category); 
                %Indices = find(L);
                
                
                if nargout >= 2
                    Blocks =  obj.Values(L);
                end
                if nargout >= 3
                    Matrix = sparse(blkdiag(Blocks{:})); %TODO: maybe introduce sparse already in Sigma.Values
                end
                if nargout >= 4
                    Vector = cell2mat(cellfun(@(x) x(:), Blocks(:), 'uniformoutput', false));
                end
            end

            
            
        end
        
        %Tab must be a cut Table of Tab_old Tab < Tab_old!!!
        function S = f_full_sigma(obj, Tab_old)
            
            Matrix_Index_cut = obj.Tab.Matrix_Index;
            Old_List_Index = obj.Tab.Old_List_Index;
            matrix_index = Old_List_Index(Matrix_Index_cut);
            
            %matrix_index = Tab_old.Matrix_Index(obj.Tab.Old_Indices);
            S = Sigma([], Tab_old.f_block(sparse(matrix_index(:,1), matrix_index(:,2), obj.f_vector, Tab_old.Numel_Datasets, Tab_old.Numel_Datasets)), Tab_old);
        end


        
        % NOT USED ANY MORE since Values is already in block structure!
        function [S, block_table] = f_block_cell(obj, N_part_inp, Spin_inp, Ener_inp, order)
            % f_block_cell(obj, N_part_inp, Spin_inp, Ener_inp, order)
            % NOT USED ANY MORE since Values is already in block structure!
            %   block_cell - method - filters sigma and groups it in blocks
            %   block_cell(N_part_inp, Spin_inp, Ener_inp, order)
            %   N_part_inp, Spin_inp and Ener_inp used to filter, can be vectors, an empty vector [] considers all 
            %   order can be set to group sigma in blocks, also defines the order, resulting blocks are returned in
            %   a cell vector S
            %   block_table illustrates the selected blocks, the selected states
            
            
            % first create a table of relevant parameters
            % define Spalten by input
            blockstructure = true;
            if   nargin < 5 || isempty(order), blockstructure = false; order = []; end
            
            
            
            
            
            %[use_N,use_E,use_S] = deal(false);
            if blockstructure
                
                if nargin < 4 || isempty(Ener_inp), Ener_inp = []; end
                if nargin < 3 || isempty(Spin_inp), Spin_inp = []; end
                if nargin < 2 || isempty(N_part_inp), N_part_inp = []; end
                [block_table,L] = obj.Energy.f_block_table( N_part_inp, Spin_inp, Ener_inp, order);
                tot_block = size(block_table,1);
                S = cell(1,tot_block);
                for n_block = 1:tot_block
                    S{n_block} = full(obj.Values(L{n_block},L{n_block}));
                end
                
                
            else
                if nargin < 4 || isempty(Ener_inp), L_E = true;
                else, L_E = ismembertol(obj.Energy.Energies, Ener_inp, obj.Tolerance); end
                if nargin < 3 || isempty(Spin_inp), L_S = true;
                else, L_S = ismember(obj.Energy.Spin, Spin_inp); end
                if nargin < 2 || isempty(N_part_inp), L_N = true;
                else, L_N = ismember(obj.Energy.N_part, N_part_inp); end
                
                
                % create logical vector:
                L = L_N & L_S & L_E & true(numel(obj.Energy.Index),1);
                
                S = obj.Values(L,L);
                block_table = L;
                
            end
            
          
        end
        
        
        
        

     
        
        function [] = f_check_entanglement(obj)
            
        end
        
        
        function [obj, analyse] = f_validate(obj, A, tolerance)
            % [obj, analyse] = f_validate(obj, A, tolerance)
            % checks if A (or, if not defined, obj.Values) fulfills the density operator properties, result is stored in analyse
            % -> calculates obj.Eigenvalues and obj.Eigenvectors
            % 
            
            if nargin < 3 || isempty(tolerance), tolerance = obj.Tolerance; end
            if nargin < 2 || isempty(A), A = obj.Values;
            else
                if ~iscell(A), A = obj.Tab.f_block(A); end
            end
            %NOTE: obj.Values is block structure, A can be matrix, vector or block_structure?
            %will be transferred into block - Tab necessary - first check structure!
            
            % analysis the following things:
            %   - checks if it is a valid density matrix (trace, positivity, hermiticity)
            %       positivity takes time!!!
            %   it calculates the eigenvalues and eigenvectors if block is smaller than 200 (Full ED)
            
            
            % check TRACE 1
            %initialize order of children:
            analyse.positivity = []; 
            analyse.is_herm = [];
            analyse.trace_one = [];
            analyse.trace = [];
            analyse.trace_positive = [];
            analyse.sum_neg_diag = [];
            analyse.sum_neg_ev = [];
            analyse.hermiticity_diff = [];
            
            
            
            analyse.trace = sum(cell2mat(cellfun(@(x) trace(x), A, 'uniformoutput', false)));
            
            if abs(analyse.trace - 1) < tolerance, analyse.trace_one = true; else, analyse.trace_one = false; end
            if all(cellfun(@(x) all(diag(full(x))) > -tolerance, A, 'uniformoutput', true)), analyse.trace_positive = true; else, analyse.trace_positive = false; end
            analyse.sum_neg_diag = sum(cellfun(@(x) sum(abs(diag(full(x))) - diag(full(x)))/2, A));
            
            % continue to work in blocks?
            
            analyse.hermiticity_diff = 0;
            analyse.sum_neg_ev = 0;
            
            if ismember('E', obj.Tab.Order)
                energy_tag = 'E'; 
                block_index = obj.Tab.Block_Index(energy_tag, {-inf, obj.Energy_Cut});
            elseif ismember('A', obj.Tab.Order)
                energy_tag = 'A';
                block_index = obj.Tab.Block_Index(energy_tag, {-inf, obj.Energy_Cut});
            else
                block_index = obj.Tab.Block_Index;
            end
            
            
            num_of_blocks = numel(block_index);
            
            
            obj.Eigenvectors = cell(num_of_blocks,1);
            obj.Eigenvalues = cell(num_of_blocks,1);
            for n_block_ind = 1:num_of_blocks
                n_block = block_index(n_block_ind);
                %if block{n_block,2} >= 1
                block = obj.Values{n_block};
                
                
                % check HERMICITY
                Herm =  block - block';
                 %ishermitian(block)
                analyse.hermiticity_diff = analyse.hermiticity_diff + sum(abs(Herm(:)));
                
                
                
                % check POSITIVITY and save eigenvalues according to NSE
                if size(block,1) < 200
                    if any(isnan(block(:))) || any(isinf(block(:)))
                        warning('inf or nan in sigma!!!')
                    else
                        [ev,ew] = eig(full(block));
                        analyse.sum_neg_ev = analyse.sum_neg_ev + sum(ew(ew < 0));
                        
                        
                        obj.Eigenvalues{n_block} = diag(ew);
                        obj.Eigenvectors{n_block} = ev;
                    end
                else
                    if ~isdeployed()
                        warning('Blocksize larger that 200 - full ed not performed')
                    end
                end
                
            end
            if  analyse.sum_neg_ev < -tolerance, analyse.positivity = false; else, analyse.positivity = true; end 
            if  analyse.hermiticity_diff < tolerance, analyse.is_herm =  true; else, analyse.is_herm = false; end
               
            obj.Analysis = analyse;
            
            %block = cell2table(block, 'variablenames', {'First_index', 'max_block_size','Index', 'Is_N_part_block', 'N_part', 'Is_Spin_block', 'Spin', 'Is_Energy_block', 'Energ'});
            
        end
        
        function [block, structure, p2, block_column_label] = f_analyse_block(obj, A, Ener_inp)
            %[block, structure, p2] = f_analyse_block(obj, A, Ener_inp)
            %
            % analyse  - method - analyses the structure of the density matrix (if there is a block structure)
            % examine blockstructure -
            % TODO: check computation time of whole function, especially of Block analysis
            % find block structure
            % problem if blocks have high periodicity (not very likely) - better algorithm?
            % for this kind of tasks sufficient - for bigger matrices implement better algorithm!!
            if nargin < 2 || isempty(A), A = obj.Values; end
            if iscell(A) || (size(A,1)~= size(A,2)), error('Blockstructure / Vectorstructure already realised (input has block / vector structure)'); end
            
            if nnz(A) < numel(A)*0.95    %check if there are at least 90 % zeros, otherwise it's not reasonable
                % maybe amd - approximate minimum degree permutation is helpful
                % have a look at amd symrcm
                
                % reordering best with symamd! - returns permutation
                % p = symamd(A)
                % return to previous ordering
                % [~,pinv] = sort(p);
                
              
               % (herm, trace and pos
                % examination may be faster)
                
                % for testing:
                %pp = randperm(size(obj.Values,1));
                %obj.Values = obj.Values(pp,pp);
                [row, column] = find(A);
                block = [];
                structure = [];    
                
                kk = 1; 
                while numel(row) >= 1
                    x = 0; y = row(1); k = 0;
                    break_flag = false;
                    while break_flag == false
                        
                        x  = column(ismember(row,y));
                        y =  unique([row(ismember(row,x)); y]);
                        if all(ismember(column,x) == ismember(row,y)), break_flag = true;  end
                        k = k+1;
                        if k > 10000, error('method inefficient / periodicity higher than 10000'); end
                    end
                    block{kk,1} =  min(y);
                    block{kk,2} = numel(y);
                    block{kk,3} = y;
              
                    % block columns 4 to 9 explains for each found sigma block whether block structure is connected to N_part, Spin or Energy
                    % questions if all states (of eigenstate basis) of a block have same value, and if individual
                    % values of all basis states are well defined (use space basis of single sector)
                    block{kk,4} = var(Ener_inp.N_part(y)) == 0 &  sum(Ener_inp.N_part_var(y)) == 0;
                    if block{kk,4}, block{kk,5} = mean(Ener_inp.N_part(y)); else, block{kk,5} = inf; end
                    block{kk,6} = var(Ener_inp.Spin(y)) == 0 & sum(Ener_inp.Spin_var(y)) == 0;
                    if block{kk,6}, block{kk,7} = mean(Ener_inp.Spin(y)); else, block{kk,7} = inf; end
                    block{kk,8} = var(Ener_inp.Energies(y)) < obj.Tolerance;
                    if block{kk,8}, block{kk,9} = mean(Ener_inp.Energies(y)) ; else, block{kk,9} = inf; end
                    
                    kk = kk+1;
                    
                    row(ismember(row,y)) = [];
                    column(ismember(column,x)) = [];
                end
                                    
                % permuation sequence, sort of blocksizes before possible
                p2 = cell2mat(block(:,3));
                %spy(obj.Values(p2,p2));
                
                
                %could be alternatively set
                % eigen_en_particle_block = sum(Ener_inp.N_part_var) == 0;
                % eigen_en_spin_block =     sum(Ener_inp.Spin_var) == 0;
                % TODO: at the moment block structure can just be true if energy eigenvalues are related to distinct particle / spin
                % values - maybe broaden this aspect

                
                if all([block{:,4}]), structure = [structure, 'N']; end
                if all([block{:,6}]) && Ener_inp.Basis.N_spin == 2, structure = [structure, 'S']; end
                if all([block{:,8}]), structure = [structure, 'E']; end
               

            else
                block = {1,size(obj.Values,1), 1:size(obj.Values,1), false, inf, false, inf, false, inf};
            end
            block_column_label = {'first_row_of_block', 'block_size', 'row_indices', 'is_block_in_N_part', 'N_part_value', 'is_block_in_Spin', 'Spin_value', 'is_block_in_Energy', 'Energy_value'};            
        end
        
        
        function [GR_sys, GK_sys] = f_greens_function(obj, Ener_inp, c_site_1, c_spin_1, c_site_2, c_spin_2, Grid, i0p, sig_tolerance)
           
            BLOCKSIZE = 5001;
            
            %NOTE: employ a counter how much weight in Green function calculation has been dropped due to energy cut...
            
            
            %function that calculates the system single particle Green's function in frequency space.
            %It uses the density matrix. It's not zero temperature and not in equilibrium.
            %there are two ways of generating a function:
            % either using a grid, or creating a function handle
            % the relevant operators also have to be chosen (ij, \sigma \sigma')
            % this function works for all combinations, but shouldn't be executed directly, maybe there is a
            % wrapper function / no parallel implementation considered so far! (just c_site_1 * c_site_2 for loops)
            En = Ener_inp;
            %TODO: adopt, big construction!
            % DONE: work in Sigma block, first term is trace( \rho(block) * Q(indices).*\frac{1}{w-\delta E} * Q^\dag(indices))
            % second block: trace( \rho(block) * Q^\dag (indices) * Q(indices) .* \frac{1}{w-\Delta E})
            % TODO: check if En.f_c delivers all necessary Q-matrices, try some configurations, make at least a
            % list!
            GR_sys = zeros(Grid.N,1);
            GK_sys = zeros(Grid.N,1);
           
            
                        
            %DONE: if sigma is blockdiagonal also in Energy / Average Energy - use those blocks and only those
            %necessary!!!
            
            
            %TODO: spinless, spin_preserving and spin_symmetric case
           %      N_part_inp  |  Spin_inp  |  c_spin
            % ___________________________________________
            % 1a)|     #       |     [ ]    |    [ ]
            % 1b)|    [ ]      |     [ ]    |    [ ]
            % 2a)|     #       |      #     |     #
            % 2b)|     #       |     [ ]    |     #
            % 2c)|    [ ]      |      #     |     #
            % 2d)|    [ ]      |     [ ]    |     #
            
            % Define modes
            if Ener_inp.Basis.N_spin == 1
                if ~isempty(c_spin_1) || ~isempty(c_spin_2)
                    warning('In spinless case, Spin_inp and c_spin should not be defined but empty!')
                end
                if ~ismember('N', Ener_inp.Structure)
                    %mode = 2;
                    %first_filter = '';
                    if ~isempty(N_part_inp)
                        warning('Case 1b, particle number not conserved, N_part_inp should not be defined but empty!')
                    end
                else
                    %mode = 1;
                    %first_filter = 'N';
                end
            else
                % with spin
                
                if ismember('S', Ener_inp.Structure)
                    %checks if right particle/spin sector is well defined
                    
                    
                    if ismember('N', Ener_inp.Structure)
                        
                        %mode = 3;
                        %first_filter = 'NS';
                       
                    else
                        %mode = 5;
                        %first_filter = 'S';
                 
                    end
                else
                    if ismember('N', Ener_inp.Structure)
                        %mode = 4;
                        %first_filter = 'N';
                     
                    else
                        %mode = 6;
                        %first_filter = '';
                     
                    end
                end
            end
            
            
          
            %NOTE: check if there is a problem, if spin_tag is zero, 
            % finding of energ_index has one argument too much
           
            
            %NOTE: table_tag shall be used for ordering, not for access of energies (take original ones!)
            
            if ismember('E', obj.Tab.Order)
                energy_tag = 'E'; 
            elseif ismember('A', obj.Tab.Order)
                energy_tag = 'A';
            else
                energy_tag = '!'; % ! is not a shortcut - will be ignored
            end
            
            %Qtab = Table(obj.Energy, ['N',spin_used, energy_block_used]);
            
            %NOTE: Sigma could be in any structure, it might be a block structure. NE, NS is possible
                      
            %NOTE: NEw - take just blocks into account which have a significant large sigma
            weight = cellfun(@(x) sum(abs(full(x(:)))), obj.Values);
            rho_block_index = obj.Tab.Block_Index(weight > sig_tolerance); 
            
            %NOTE: obj.Tab stores the block structure of Sigma! surf through all relevant blocks.
            
            %DONE: Take only those energies into account where sigma is significant!!! (sig_tolerance used for diagonal)
            %DONE: check if energy_cuts are stored in sigma object
           
            %disp([num2str(numel(energ_index)), ' Blocks are used in order to calculate the retarded Green function'])
           
            %TODO: make ordered for loop over all relevant entries in sigma, order by weight.
            % Check if weight is member of obj.Tab.Categ
            
            
            
            for ind_m = 1:numel(rho_block_index)
                n_block = rho_block_index(ind_m);
                rho = full(obj.Values{n_block});
                [rho_m, rho_n] = size(rho);
                    
                %TODO: Question: are we checking if Sigma is blockdiagonal in N and S or Energy?
                if ismember('N', obj.Tab.Order)
                    n_part = obj.Tab.N_part(n_block);
                    %L_N = obj.Tab.N_part == n_part;
                else
                    n_part = [];
                    %L_N = true;
                end
                
                
                if ismember('S', obj.Tab.Order)
                    spin = obj.Tab.Spin(n_block);
                    %L_S = obj.Tab.Spin == spin;
                else
                    spin = [];
                    %L_S = true;
                end
                
                %n_block_first = find(L_N & L_S, 1, 'first'); %TODO: check if this works
                %block_index = sum(obj.Tab.Block_Size(n_block_first:(n_block-1))) + (1:obj.Tab.Block_Size(n_block));
                
                %TODO: now find out according to which indexing scheme we are proposing, how is Tab.Order
                %related to Ener_inp.Structure. Define relations as Q_mat is derived from Ener_inp and Sigma
                %structure could be more general (non secular approximation)
                
                %Draw two blocks, one energy and introduce Average Energy
                
                
                
                %NOTE: The list indices (those indices in Energy class corresponding to a block in the Table)
                % are returned by the List_Index function
                % to know which parts of the Qmatrices one should take, we have the SubIndex function to search
                list_index = obj.Tab.Old_List_Index('I', n_block);
                
                
                E_n = En.Energies('I',list_index);
                sub_block_index = En.Subindex('NS',  {n_part, spin}, 'I', {list_index});
                %E_m = En.g_En(n_part,spin);
                %E_m = E_m(block_index);
                %E_n = E_m;
                
                %if rho > 0
                %    disp(rho)
                %    disp(ind_m)
                %end
                
                %for ind_n = 1:numel(obj.Tab.Block_Index)
                dag = 1;
                %distinguish between cases, how can Tab be organized?
                % in order to evaluate single particle Green's function we can shorten our
                % calculation by working in distincive particle spin sectors since c^\dag c just work
                %1) NSE
                %2) NE
                %3) N
                % according to c dagger operator work the following:
                % <      c^s_{c_spin, c_site} | N, S >
                
                %[E_p, E_q] = deal(0);
                % NOTE: preinitialized with zero if some part in c_function is not applied
                % TODO: Use flags to indicate C function which parts shall be calculated to prevent this
                C1 = zeros(rho_n,1);
                C2 = zeros(1,rho_m);
                C3 = zeros(rho_n,1);
                C4 = zeros(1,rho_m);
                E_p = 0;
                E_q = 0;
                creation = false;
                annihilation = false;
                %rho = full(obj.g_matrix('NS', n_part, spin));
                
                % NOTE: creation part:  rho \hat(c_i) c_j^dag
                if (n_part + dag >= 0) && (n_part + dag <= Ener_inp.Basis.N_site*Ener_inp.Basis.N_spin) && ...
                        ( ~ismember('S', Ener_inp.Structure) ||   ismember(spin + dag * c_spin_2, Ener_inp.Basis.f_spin_range(n_part + dag)) )
                    
                    %E_p = En.g_En(n_part+dag,spin + dag*c_spin_2);
                    %NOTE: Access of energies: take original ones (not averaged!)
                    %NOTE: not defined conserved quantities will be [] and ignored since [] + 1 = []
                    E_p = En.Energies(['NS', energy_tag], n_part + dag, spin + dag*c_spin_2, {-inf, obj.Energy_Cut_2});
                    if ~isempty(E_p)
                        % right state first, left state second
                        
                        creation = true;
                        %NOTE: this works since not defined conserved quantities will be ignored since they are []
                        energy_subindex = En.Subindex('NS', {n_part + dag, spin + dag*c_spin_2}, energy_tag, {{-inf, obj.Energy_Cut_2}});
                        %energy_subindex = En.f_energy_subindex(n_part + dag, spin + dag*c_spin_2, {-inf, obj.Energy_Cut_2});
                        
                        %DONE: use g_Qmat with respecting energy cuts!!!
                        C2 = En.g_Qmat(dag,n_part,spin,c_site_2,c_spin_2,obj.Energy_Cut_2); %c_j^dag
                        C2 = full(C2(energy_subindex, sub_block_index));
                        C1 = En.g_Qmat(-dag, n_part + dag, spin + dag * c_spin_2, c_site_1, c_spin_1,obj.Energy_Cut_2); %c_i
                        C1 = full(C1(sub_block_index, energy_subindex));
                        
                    else
                        E_p = 0;
                        warning('creation part empty / out of energy cut')
                        
                    end
                end
                %NOTE: annihilation part: rho c_j^dag \hat(c_i)
                if (n_part - dag >= 0) && (n_part - dag <= Ener_inp.Basis.N_site*Ener_inp.Basis.N_spin) && ...
                        ( ~ismember('S', Ener_inp.Structure) || ismember(spin - dag * c_spin_1, Ener_inp.Basis.f_spin_range(n_part - dag)) )
                    % right state first, left state second, N-dag sector
                    
                    %NOTE: not defined conserved quantities will be [] and ignored since [] + 1 = []
                    E_q = En.Energies(['NS', energy_tag], n_part - dag, spin - dag*c_spin_1, {-inf, obj.Energy_Cut_2});
                    %E_q = En.g_En(n_part - dag, spin - dag*c_spin_1);
                    if ~isempty(E_q)
                        
                        annihilation = true;
                        %energy_subindex_2 = En.f_energy_subindex(n_part - dag, spin - dag*c_spin_1, {-inf, obj.Energy_Cut_2});
                        
                        %NOTE: this works since not defined conserved quantities will be ignored since they are []
                        energy_subindex_2 = En.Subindex('NS', {n_part - dag, spin - dag*c_spin_1}, energy_tag, {{-inf, obj.Energy_Cut_2}});
                        
                        C4 = En.g_Qmat(-dag, n_part, spin, c_site_1, c_spin_1, obj.Energy_Cut_2); %c_i
                        C4 = full(C4(energy_subindex_2, sub_block_index));
                        C3 = En.g_Qmat(dag, n_part - dag, spin - dag*c_spin_1, c_site_2, c_spin_2, obj.Energy_Cut_2); %c_j^dag
                        C3 = full(C3(sub_block_index, energy_subindex_2));
                        
                    else
                        E_q = 0;
                        warning('annihilation part empty / out of energy cut')
                        
                    end
                end
                %disp(ind_m)
                if annihilation || creation
                    
                    
                    
                    points = Grid.Points;
                    %disp(['energy_block: ', num2str(ind_m), ', creation: ', num2str(creation), ', annihilation: ', num2str(annihilation)])
                    if isempty(C1) || isempty(C2) || isempty(C3) || isempty(C4)
                        error('C function cannot be called with empty array - something went wrong!!!')
                    end
                    
                    
                    % NOTE: generated with Matlab Coder from m_GF.m (the lines below), as fast as my own code,
                    % double as fast as compared to Matlab code
                    if isunix
                        [GR, GK] = m_GF_mex(points, i0p, E_p, E_n, E_q, complex(C1), complex(C2), complex(C3), complex(C4), complex(rho));
                    else
                        [GR, GK] = m_GF(points, i0p, E_p, E_n, E_q, complex(C1), complex(C2), complex(C3), complex(C4), complex(rho));
                    end
                    % NOTE: my own C - code (does not work well - some errors - could not figure out why...)
                    %tic
                    %[GR, GK] = c_GF(points, i0p, E_p, E_n, E_q, complex(C1), complex(C2), complex(C3), complex(C4), complex(rho));
                    %toc
                    
                    % NOTE: easy simple implementation in Matlab (3D array operations)
                    %Test_R = diagsum(mmat(mmat(rho, bsxfun(@times, (1./(permute(points,[2,3,1]) - (E_p' - E_n) + i0p*1i)),C1)), C2) + ...
                    %   mmat(mmat(rho, C3), bsxfun(@times, (1./(permute(points, [2,3,1]) - (E_n' - E_q) + i0p*1i)), C4)),1,2);
                    
                    %Test_G = diagsum(-2i*i0p*   mmat(rho, mmat(bsxfun(@times, (1./((permute(points,[2,3,1]) - (E_p' - E_n)).^2 + i0p^2)), C1), C2)) + ...
                    %                 2i*i0p*   mmat(rho, mmat(C3, bsxfun(@times, (1./((permute(points,[2,3,1]) - (E_n' - E_q)).^2 + i0p^2)), C4))), 1,2);
                    
                    % NOTE: m_GF.m: test case
                    %tic
                    %[E_pp, E_nn] = meshgrid(E_p, E_n);
                    %[E_mm,E_qq] = meshgrid(E_n, E_q);
                    %GR = complex(zeros(numel(points),1));
                    %GK = complex(zeros(numel(points),1));
                    %k = find((abs(points + 1.486)) < 10^-4)
                    %for k = 1:numel(points)
                    %   GR(k) = trace( rho * (1./(points(k) - (E_pp - E_nn)+i0p*1i).*C1) * C2 + rho * C3 * (1./(points(k) - (E_mm - E_qq) + i0p*1i) .* C4));
                    %   GK(k) = trace( -2i*i0p * rho * (1./( (points(k) - (E_pp - E_nn)).^2 + i0p^2) .* C1) * C2  + 2i*i0p * rho * C3 * (1./( (points(k) - (E_mm - E_qq)).^2 + i0p^2) .* C4) );
                    %end
                    %toc
                    
                    %figure
                    %subplot(2,2,1), plot(points, real(GK)), hold on, plot(points, real(Test_G))
                    %subplot(2,2,3), plot(points, imag(GK)), hold on, plot(points, imag(Test_G))
                    %subplot(2,2,2), plot(points, real(GR)), hold on, plot(points, real(Test_R))
                    %subplot(2,2,4), plot(points, imag(GR)), hold on, plot(points, imag(Test_R))
                    
                    
                    %if any(imag(GR) > 0)
                    %    warning('should be smaller than 0')
                    %end
                    GR_sys = GR_sys + GR;
                    GK_sys = GK_sys + GK;
                end
                
            end
            
            
            
            
            
            %DONE: Think about Green function class and inplement also a Keldysh calculation!
            %TODO: Check if parts of calculation can be recycled (faster) or if an extra function is wise
            
            
        end
        
        function [GR,GK] = f_greens_function_full(obj, Ener_inp, grid, spin_symmetric, i0p)
            % TODO: spin_symmetric: indicates if Green functions are spin_symmetric (check in future bath class version
            if nargin < 4 || isempty(spin_symmetric), spin_symmetric = false; end
            if nargin < 5 || isempty(i0p), i0p = grid.Imag_Infinitesimal/10; end %TODO: make it bigger?
            grid_N = grid.N;
            %NOTE: check which kind of operators in the central sysmtem have to be formed:
            N_site = Ener_inp.Basis.N_site;
            N_spin = Ener_inp.Basis.N_spin;
            
            sig_tolerance = 10^-6;
            
            disp(['Creation of Green function, weight of blocks considered: ', num2str(sum(obj.f_diagonal_weight(abs(obj.f_diagonal_weight)> sig_tolerance))), ',   used tolerance: ', num2str(sig_tolerance)])
            
            
            k = 0;
            %DONE: implement for Spin = 1
            %DONE: implement spin symmetric case!!!
            
            % Think about cases:
            % 1) spinless
            % Structure of GR and GK:
            %   (N_site, N_site, omega_points)
            % 2) with spin
            %   (2*N_site, 2*N_site, omega_points)
            % 2a) spin conserved:
            %    just combinations with same spin
            %    2aa) spin symmetric
            %
            %    2ab) not spin symmetric
            %
            % 2b) spin not conserved:
            %    every variation of spin_up and spin_down
            %    2ba) spin symmetric
            %
            %    2bb) not spin symmetric
            
            
            
            
            %NOTE: if Hamiltonian preserves the spin then < c_up c^dag_down > is zero!!!
           
            
            
            %TODO: implement spin dependence of coupling and bath!!!
            if N_spin == 1
                spin_1_id_range = 1;
                spin_2_id_range = 1;
                copy_spin = false;
            else
                %N_spin == 2
                spin_1_id_range = 1:2;
                spin_2_id_range = 1:2;
                
                if spin_symmetric
                    copy_spin = true;
                    %NOTE: spin_symmetric is just implemented for spin conserved Hamiltonians (no c_1up c_1down entries. 
                    % Thus don't loop through spin_1 = 1 (or spin_2 = 1), this will be done later
                    spin_1_id_range = 1;
                    
                    if ~ismember('S', obj.Block_structure)
                        warning('spin symmetric for not spin conserved objects not implemented!')
                    end
                else
                    copy_spin = false;
                end
            end
            
            [GR,GK] = deal(zeros(N_spin * N_site, N_spin * N_site, grid_N));

            for c1 = 1:N_site
                for spin_1_id = spin_1_id_range
                    if ismember('S',obj.Block_structure)
                        %spin conserved, there are just entries with the same spin (G is block diagonal in spin)
                        spin_2_id_range = spin_1_id;
                    end
                    for c2 = 1:N_site
                        for spin_2_id = spin_2_id_range
                            k = k+1;
                            if k == 1, check_time = tic; end
                            %disp(num2str(k))
                            if N_spin == 1
                                spin_1 = [];
                                spin_2 = [];
                            else
                                spin_1 = (spin_1_id-1.5)*2;
                                spin_2 = (spin_2_id-1.5)*2;
                            end
                       
                            [g_r,g_k] = obj.f_greens_function(Ener_inp, c1, spin_1, c2, spin_2, grid, i0p, sig_tolerance);
                            
                            if ~isdeployed()
                                if k == 1, time_pre = toc(check_time); 
                                    disp(['Estimation for calculation of full Greens function: ', num2str(time_pre*(N_site^2*numel(spin_1_id_range) * numel(spin_2_id_range)))]); 
                                else
                                %    disp([num2str(100*k/(N_site^2*numel(spin_1_range) * numel(spin_2_range))), '% calculated'])
                                end
                            end
                            GR(c1 + N_site*(spin_1_id-1), c2 + N_site*(spin_2_id-1),:) = g_r;
                            GK(c1 + N_site*(spin_1_id-1), c2 + N_site*(spin_2_id-1),:) = g_k;
                        end
                    end
                end
            end
            
            %NOTE: copy spin to spin up area since Hamiltonian is spin preserving ans spin symmetric!
            if copy_spin
                GR((1:N_site) + N_site, (1:N_site) + N_site, :) = GR(1:N_site, 1:N_site,:);
                GK((1:N_site) + N_site, (1:N_site) + N_site, :) = GK(1:N_site, 1:N_site,:);
            end
            
            
        end
        
        
        function [vector] = f_vector(obj, varargin)
            if nargin == 1
                vector = obj.Tab.f_vector(obj.Values);
            else
                [~,~,~, vector] = obj.g_superblocks(varargin{:});
            end
        end
        
        function [matrix] = f_matrix(obj, varargin)
            if nargin == 1
                matrix = obj.Tab.f_matrix(obj.Values);
            else
                [~,~,matrix] = obj.g_superblocks(varargin{:});
            end
        end
        
        function block = f_block(obj, varargin)
            if nargin == 1
                block = obj.Tab.f_block(obj.Values);
            else
                [~, block] = obj.g_superblocks(varargin{:});
            end
        end
        
        function weight = f_diagonal_weight(obj,L)
            if nargin < 2 || isempty(L)
                weight = cellfun(@(x) sum(abs(diag(full(x)))), obj.Values);
            else
                weight = cellfun(@(x) sum(abs(diag(full(x)))), obj.Values);
                weight = weight(L);
            end
        end
            
        function diagonal = f_diagonal(obj)
            diagonal = diag(full(obj.g_matrix));
        end
        
        % NOT USED at the moment
        function [S] = Ff_block(obj, N_part_inp, Spin_inp, Ener_inp)
            % f_block(obj, N_part_inp, Spin_inp, Ener_inp)
            % NOT USED at the moment
            % returns the block to particle number and/or spin sector and/or Energy
            
%             if nargin < 3 || isempty(Ener_inp), L_E = true; 
%             else L_E = abs(Ener_inp - obj.Energy.En) < obj.Tolerance; end
%             if nargin < 2 || isempty(Spin_inp), L_S = true; 
%             else L_S = Spin_inp == obj.Energy.Spin; end
%             if nargin < 1 || isempty(N_part_inp), L_N = true; 
%             else L_N = N_part_inp == obj.Energy.N_part; end
%             
            %if vector style is allowed
            if nargin < 4 || isempty(Ener_inp), L_E = true; 
            else, L_E = ismembertol(obj.Energy.Energies, Ener_inp, obj.Tolerance); end
            if nargin < 3 || isempty(Spin_inp), L_S = true; 
            else, L_S = ismember(obj.Energy.Spin, Spin_inp); end
            if nargin < 2 || isempty(N_part_inp), L_N = true; 
            else, L_N = ismember(obj.Energy.N_part, N_part_inp); end
            
            
            % create logical vector:
            L = L_N & L_S & L_E ;
            
            S = obj.Values(L,L);
        end
        
        
        function percent = f_exp_val(obj, En, Tab, code)
            %NOTE: valid in the case of a spin and particle preserving Hamiltonian!
            
            %NOTE: code: [nan, 0, 1, nan, nan, nan] == < c_1u^dag c_2u >
            % (as defined in Bas.Binary_Label)
            
            
%             dom = obj.s_dom_blocks;
%             
%              %  dom.Old_Indices
%             
%             dom.weight;
%             
             binar = En.Basis.f_bin;
            L_full = code == 1;
            L_empty = code == 0;
            
            old_index = obj.Tab.Old_List_Index();
            
            %operator in the single space basis
            L = all(binar(:, L_full) == 1,2) & all(binar(:, L_empty) == 0,2);
            
            %NOTE: Table to manage Basis space
            Tab_NS = Table(En, 'NS');
            %NOTE: returns all relevant particle spin sectors - may be adopted to other conservation laws!
            tab = unique([obj.Tab.N_part, obj.Tab.Spin], 'rows');
            
            %which particle and spin sectors of the space basis are used?
            L_relevant = cell2mat(arrayfun(@(x,y) Tab_NS.List_Index('NS', x,y), tab(:,1), tab(:,2), 'uniformoutput', false));
            
            %block_size vector of all relevant blocks
            block_size_NS = arrayfun(@(x,y) Tab_NS.Block_Size('NS', x,y), tab(:,1), tab(:,2));
            %shrinkage to relevant blocks and packaging in a cell
            L_NS = L(L_relevant);
            L_cell = mat2cell(L_NS, block_size_NS,1);
            
            
            %eigenvector to transform rho to space basis ev * rho * ev'
            %rows correspond to space basis (NS), columns to many body basis (NSI)
            EV = arrayfun(@(x,y) En.Eigenvector(En.Index('NS',x,y), En.Index('INS', old_index, x,y)),  tab(:,1), tab(:,2), 'uniformoutput', false) ;
            
            %map rho to NS structure to match eigenvectors (INS) - I through cut table
            rho_cell = arrayfun(@(x,y) obj.f_matrix('NS', x,y),  tab(:,1), tab(:,2), 'uniformoutput', false) ;
            
            %calculate value
            percent = sum(cellfun(@(x,y,z) trace(x * z * x' * diag(y) ), EV, L_cell, rho_cell));
            percent = denoise(percent, 10^-14);
            
            %% Test via full matrix operation in full space - to check
%             T = zeros(En.Numel_En, En.Numel_En);
% 
%             for k = 1:numel(obj.Values)
%                 T(obj.Tab.Li_Index{k},obj.Tab.Li_Index{k}) =  obj.Values{k};
%             end
%             ev = En.Eigenvector;
%             
%             percent = trace( ev*T * ev' * diag(double(L)));
  
            
        end
            
        function full_current_map = f_full_current(obj, En)
            vec = nan(1, En.Basis.N_site * En.Basis.N_spin);
            full_current_map = zeros(length(vec));
            for k = 1:(En.Basis.N_site * En.Basis.N_spin)
                for j = 1:(En.Basis.N_site * En.Basis.N_spin)
                    t = vec;
                    t(k) = 1; 
                    t(j) = 0;
                    full_current_map(k,j) = obj.f_exp_val(En,obj.Tab, t);
                end
                
                
                
            end
        end
        
        
        function index = f_index(obj,varargin)
            
            index = obj.Tab.Block_Index(varargin{:});
            
        end
        
        
%             % f_vector(obj, N_part_inp, Ener_inp, Spin_input)
%             % TODO: handle options, define what this function shall be able to do
%             % at the moment capable of transforming into vector (same as with Table.f_vector (according to it's structure / filter) 
%             % for application to superoperators
%             blocktable = obj.Energy.f_block_table([],[],[],'NS');
%             blocks = obj.f_block_cell([],[],[],'NS');
%             block_temp = cellfun(@(x) x(:),blocks,'uniformoutput',false);
%             block_vec = cell2mat(block_temp(:));
%            
%             disp([(1:numel(block_vec))',block_vec])
%             %36 Zeilen...
%             %herausfinden, welche davon für bestimmtes $\kappa s$ relevant sind
%             % -> mithilfe von Qtable
%             QQ = obj.Energy.f_create_Q_table;
%            % zB dagger = +1, spin = -1 
%            dagger = 1; 
%            spin = -1;
%             table = QQ.s_print(dagger, spin);
%             L1 = table(:,4) ~= 0;
%             sig_ind = table(:,3).^2;
%             lower = [0; cumsum(sig_ind(1:end-1))]+1; 
%             upper = sig_ind-1 + lower;
%             
%             sig_index = arrayfun(@(x,y,z) repmat((x:y)',z,1),lower(L1), upper(L1), table(L1,4),  'UniformOutput', false)
%             numel(cell2mat(sig_index(1:end)))
%             sig_index{:}
%             
%         end
        %TODO: write new
        function  [] = s_plot(obj, varargin)
            
            
            % selecting
            % a) by energy     b) by weight - done via varargin input
            block_index_orig = obj.Tab.Block_Index(varargin{:});
            if ~ismember('W', obj.Tab.Categ.f_shortnames)
                weight = cellfun(@(x) sum(abs(diag(full(x)))), obj.Values);
                obj.Tab = obj.Tab.f_add_categ({'Weight', 'W'}, weight);
            end
            if ~ismember('X', obj.Tab.Categ.f_shortnames)
                xweight = cellfun(@(x) sum(sum(abs(full(x) - diag(diag(x))))), obj.Values);
                obj.Tab = obj.Tab.f_add_categ({'XWeight', 'X'}, xweight);
            end


            fig = figure;
            plot_settings(false)
            tabgrp = uitabgroup(fig);
            sort_text = {'Tab structure', 'Energy', 'Weight'};
            if ~ismember('N', obj.Tab.Order) || ~ismember('E', obj.Tab.Order)
                k_max = 1;
            else
                k_max = 3;
            end
            
            for k = 1:k_max
                tab(k) = uitab(tabgrp, 'Title', sort_text{k}, 'backgroundcolor', 'w');
                
                % sorting - 3 tabs!
                % a) by energy
                if k == 1
                    ix = 1:numel(block_index_orig);
                elseif k == 2
                    if ismember('A', obj.Tab.Order)
                        energ = obj.Tab.Average_Energies(block_index_orig);
                    else
                        energ = obj.Tab.Energies;
                    end
                    
                    [~,ix] = sort(energ);
                elseif k ==3
                    
                    % b) by weight
                    weight = obj.Tab.Weight(block_index_orig);
                    [~, ix] = sort(weight,'descend');
                end
                
                block_index = block_index_orig(ix);
                
                
                % plotting
                % a) via f_matrix
                
                F = full(obj.f_matrix(block_index));
                
                
                % b) via List Index
                %list_index = obj.Tab.List_Index(block_index);
                %F = full(obj.f_matrix());
                %F = F(list_index, list_index);
                
                
                
                axes('parent', tab(k))
               
                
                h = guidata(gcf);
                h.button = uicontrol('style','togglebutton',...
                        'Units', 'normalized', 'position', [0.82,0.30,0.089,0.024], ...
                        'string', 'Axis equal');
                set(h.button, 'callback', {@rescale_complex_plot});
                guidata(gcf,h);
                
                plot_complex(F,k);
                %image(abs(F), 'CDataMapping','scaled')
                min_max = max(abs(F(:)));
                plot_colorbar;
                set(gca,'clim', [-min_max,min_max])
                set(gcf, 'color', 'w')
                


%                 
%                 btn = uibutton(fig,'state',...
%                'Text', 'Record',...
%                'Value', true,...
%                'Position',);
%                 % mark blocks
                hold on
                if ismember('N', obj.Tab.Order)
                    obj.Tab.s_plot(block_index)
                end
                title(['trace: ', num2str(sum(obj.Tab.Weight(block_index))), ', off diag: ', num2str(sum(obj.Tab.XWeight(block_index)))])
                % label blocks
                
               
                ax=gca();
                %% Link the 'Limits' property
                if ismember('N', obj.Tab.Order)
                   r1(k)=ax.YAxis(1);
                   r2(k)=ax.YAxis(2);
                   setappdata(ax, 'YLim_listeners', linkprop([r1(k) r2(k)],'Limits')); 
                end
                %h(k) = linkprop([r1(k) r2(k)],'Limits');
                
                
                
            end
            
            
            
        end
        
        function  [] = s_plot_paper(obj, Tab_N)
            
            set(0, 'defaultTextInterpreter', 'latex'); 
            set(0, 'DefaultAxesFontSize', 12)
            % selecting
            % a) by energy     b) by weight - done via varargin input
            block_index_orig = obj.Tab.Block_Index();
            if ~ismember('W', obj.Tab.Categ.f_shortnames)
                weight = cellfun(@(x) sum(abs(diag(full(x)))), obj.Values);
                obj.Tab = obj.Tab.f_add_categ({'Weight', 'W'}, weight);
            end
            if ~ismember('X', obj.Tab.Categ.f_shortnames)
                xweight = cellfun(@(x) sum(sum(abs(full(x) - diag(diag(x))))), obj.Values);
                obj.Tab = obj.Tab.f_add_categ({'XWeight', 'X'}, xweight);
            end


            fig = figure;
            plot_settings(false)
            
            sort_text = {'Tab structure', 'Energy', 'Weight'};
            if ~ismember('N', obj.Tab.Order) || ~ismember('E', obj.Tab.Order)
                k_max = 1;
            else
                k_max = 3;
            end
            
            
                
                
                % sorting - 3 tabs!
                % a) by energy
                
                ix = 1:numel(block_index_orig);
                
                
                block_index = block_index_orig(ix);
                
                
                % plotting
                % a) via f_matrix
                
                F = full(obj.f_matrix(block_index));
                
                
                % b) via List Index
                %list_index = obj.Tab.List_Index(block_index);
                %F = full(obj.f_matrix());
                %F = F(list_index, list_index);
                
                
                
                
               
                
                h = guidata(gcf);
                h.button = uicontrol('style','togglebutton',...
                        'Units', 'normalized', 'position', [0.82,0.30,0.089,0.024], ...
                        'string', 'Axis equal');
                set(h.button, 'callback', {@rescale_complex_plot});
                guidata(gcf,h);
                
                plot_complex(F,1);
                %image(abs(F), 'CDataMapping','scaled')
                min_max = max(abs(F(:)));
                plot_colorbar([],[],[],true);
                colorbar('location', 'south')
                set(gca,'clim', [-min_max,min_max])
                set(gcf, 'color', 'w')
                


%                 
%                 btn = uibutton(fig,'state',...
%                'Text', 'Record',...
%                'Value', true,...
%                'Position',);
%                 % mark blocks
                hold on
                if ismember('N', obj.Tab.Order)
                    Tab_N.s_plot_colorwidth('b',5)
                    obj.Tab.s_plot_colorwidth('r',3)
                end
                %title(['trace: ', num2str(sum(obj.Tab.Weight(block_index))), ', off diag: ', num2str(sum(obj.Tab.XWeight(block_index)))])
                % label blocks
                
               
                
                
                
          
            
            
            
        end
        
        
        function A = s_dom_blocks(obj, limit, blocks, En)
            % s_dom_blocks(obj, limit)
            % shows dominant blocks, considers maximum diagonal entry of each block as criteria
            % standard limit is 10^-6;
            if nargin < 2 || isempty(limit), limit = 10^-10; end
            if nargin < 3 || isempty(blocks), blocks = true; end
            if nargin < 4 || isempty(En)
                if blocks == false && isempty(obj.Tab.Energy)
                    error('in Sigma.s_dom_blocks: no Energy class available for nonblock output');
                else
                    En = obj.Tab.Energy;
                end
            end
            
            if blocks
                ev = cellfun(@(x) denoise(eig(full(x)), obj.Tolerance), obj.Values, 'uniformoutput', false);
                weight = cellfun(@(x) sum(abs(x)), ev);
                neg_ev = cellfun(@(x) sum(x(x<0)), ev); 
                L_show = abs(weight) > limit;
                weight_tab = array2table(weight, 'variablenames', {'ev_weight'});
                neg_ev_tab = array2table(neg_ev, 'variablenames', {'neg_ev'});
                
                A = [obj.Tab.s_table, weight_tab, neg_ev_tab];%(weight, 'weight');
                A = A(L_show,:);
                
                %sort
                [~, ix] = sort(weight(L_show));
                A = A(flipud(ix),:);
            else
                ev = denoise(eig(full(obj.f_matrix)), obj.Tolerance);
                weight = abs(ev);
                neg_ev = sum(ev(ev<0));
                L_show = abs(weight) > limit;
                weight_tab = array2table(weight, 'variablenames', {'ev_weight'});
                neg_ev = array2table(neg_ev, 'variablenames', ['neg_ev']);
                
                A = [En.s_print, weight_tab, neg_ev];
                [~, ix] = sort(weight(L_show));
                
                A = A(flipud(ix),:);
                
                
            end
            
            
            
        end
        
    end
    
    
    
end



