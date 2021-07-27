classdef Table
    %author: Gerhard Dorn
    %date: dezember 15
    %use: store lookup table for any kind of structure, enable indexing, and conversion of particle/spin/energy sector to
    %corresponding indices, parent class (children: Q_Table, Filter_Table)
    %features: -
    %          -
    %          -
    
    properties
        %TODO: verify if original Energy (Data class) has to be stored here
        Energy
        Order
        Categ
        
%         N_part
%         Spin
%         Energies
%         Block_Size
%         Index
        %refers to last indices - generated when table cut is used, 
        % otherwise created dynamically
        Li_Index
        
        % derived from current Table class - storage not needed - think about dynamically store entries (just one
        % calculation
        
        Vec_Index
        Diag_Index
        
        Is_Cut % True of false - indicates whether there is a cut table or not
        
        % for counting states:
        Numel_Blocks % number of blocks according to which the data is grouped
        Numel_Datasets %represents the number of rows
        Numel_Vector_Elements % represents the number of vector elements

    end
    
    
    methods
        function obj = Table(energy_inp, order_inp, Categ_Tab, vector_index_inp, diag_index_inp, n_vector_elements)
            % two constructors
            % version a) Table(energy_inp, order_inp) - creates Table according to order_inp and Energy class
            % version b) Table(... ) all parameters - if all are known, all are written (useful for subclasses)
            
            % TODO: make Table more generell - what do I need:
            % - a data set - given by energy_inp
            %   what do I use:
            %   - Structure (derived from Categ)
            %   - Categ
            %   - 
            %check input parameters
            
            %NOTE: order_inp necessary since we want to have rougher sortings (non secular approximation etc.)
            
            
            
            if nargin < 1 || isempty(energy_inp), error('Table needs Energy class to work properly'); end
            %NOTE: some standard assignment - not so meaningful at all since index is always included 
            %NOTE: we try to avoid energy.Structure since it's not a must have! - should be applicable to any
            %object containing a Category class.
            if nargin < 2  order_inp = energy_inp.Structure; %NOTE: take the structure from Energy class
                % which was defined via the Hamiltonian
            else
                if ~ischar(order_inp) || ~all(ismember(order_inp, energy_inp.Categ.f_shortnames))
                    error(['Table can be ordered according to ', energy_inp.Categ.f_shortnames]);
                end
%                 
            end
            
            if nargin < 3
                %create order
                
                % NOTE: important part for grouping!
                data_columns = cell2mat(values(energy_inp.Categ.List,num2cell(order_inp)));
                %creates a table of the categories, prior step to use unique - rows
                Ener_tab = [energy_inp.Categ.Data(data_columns).Values];
                
                %NOTE: adds a Table with one block
                if isempty(data_columns)
                    Ener_tab = ones(energy_inp.Numel_En,1); 
                end
                
                %DONE: perform rounding procedure already done in Category class - there is no sense in storing
                %values with a given error without having rounded them already
                
                %TODO:  question: - is this legimite to round already? is this more precise? or should there be
                %just a grouping? small way to quasidegenerate states...
                % are energy differences important, so it's good to keep them? what's the cost for positivity?
                
                
                
                [Block_table,~, frequ] = unique(Ener_tab,'rows');
                block_size_inp = histc(frequ,1:size(Block_table,1));
                
                %TODO: make the allocation of the new categories dynamic (take the old names and shortnames and
                %just add them.
                
                % now create the new Category
                Categ_Tab = Category();
               
                n_blocks = size(Block_table,1);
                %The index for the table:
                Categ_Tab = Categ_Tab.f_append({'Index', 'I'}, (1:n_blocks).') ;
                
                long_name = energy_inp.Categ.f_longnames(data_columns);
                %TODO: also take the Tolerance and Variance
                for k = 1:numel(order_inp)
                    Categ_Tab = Categ_Tab.f_append({long_name{k}, order_inp(k)}, Block_table(:,k));
                end
                
                
                %NOTE: add extra blocksize category
                Categ_Tab = Categ_Tab.f_append({'Block_Size', 'B'}, block_size_inp);
                
                
                
                %TODO: think about vector_index and its generalization, always derived from current table
                % independent from table cut.
                % maybe implement dynamic storage
                % pregeneration not so useful...
                
                if energy_inp.Categ.Data(1).Numel > 7000
                    if ~isdeployed()
                        warning('Creation of vector index may create a "out of memory"')
                    end
                    vector_index_inp = {};
                    diag_index_inp = {};
                    n_vector_elements = sum(block_size_inp.^2);
                else
                    deg2 = cumsum([0; block_size_inp(1:end-1).^2]);
                    vector_index_inp = arrayfun(@(x,y) (y + (1:(x^2)))',block_size_inp, deg2, 'UniformOutput', false);
                    n_vector_elements = deg2(end)+ block_size_inp(end).^2;
                    %TODO: think about diag_index and its generalization
                    diag_index_inp = arrayfun(@(x,y) (y + (1:(x+1):(x^2))'),block_size_inp, deg2, 'UniformOutput', false);
                end
                %obj = Table(energy_inp, order_inp, n_part_inp, spin_inp, energies_inp, block_size_inp, index_inp, diag_index_inp);
                if ~isempty(data_columns) && var([energy_inp.Categ.Data(data_columns).Numel]) ~= 0
                    if ~isdeployed()
                        warning('Data sets of given Categories have not the same number of elements!')
                    else
                        disp('Warning: Data sets of given Categories have not the same number of elements!')
                    end
                end
                

                
                
            else %constructor for subclass
                %discuss whether use of Vec_Index and Diag_Index are useful!!!
                
                %TODO: define what happens when Vec_Index and Diag_Index are not set!!!
                
            end
            %fill in (use also for subclasses)
            obj.Energy = energy_inp;
            obj.Order = order_inp;
            
            if isempty(data_columns)
                n_datasets = energy_inp.Numel_En;
            else
                n_datasets = mean([energy_inp.Categ.Data(data_columns).Numel]);
            end
            n_blocks = size(Categ_Tab.Data(1).Values,1);
            obj.Categ = Categ_Tab;
            obj.Numel_Datasets = n_datasets;
            obj.Numel_Blocks = n_blocks;
            
            %TODO: may create an overflow error
            obj.Numel_Vector_Elements = n_vector_elements;
            obj.Vec_Index = vector_index_inp; % cell
            obj.Diag_Index = diag_index_inp;    % cell
            
            obj.Is_Cut = false;
            obj.Li_Index = [];
                       
            % implement Index is Order and Order in Category Sequence
            % #1 is always Index (since it is always defined)
            
            
        end
        
        
        function obj = f_add_categ(obj, name, values)
            obj.Categ = obj.Categ.f_recreate_lists;
            obj.Categ =  obj.Categ.f_append(name, values);
            
        end
        
        function index = Index(obj, varargin)
            error('Use Block_Index instead, Index not used any more')
            index =  obj.Categ.f_access('Index', varargin{:});
        end
        
        function index = Block_Index(obj, varargin)
            index =  obj.Categ.f_access('Index', varargin{:});
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
        
        function block_size = Block_Size(obj, varargin)
            block_size = obj.Categ.f_access('Block_Size', varargin{:});
        end
        
        function energ = Average_Energies(obj, varargin)
            energ =  obj.Categ.f_access('Average_Energies', varargin{:});
        end
        
        function old_indices = Old_Indices(obj, varargin)
            old_indices =  obj.Categ.f_access('Old_Indices', varargin{:});
        end
        
        function weight = Weight(obj, varargin)
            weight =  obj.Categ.f_access('Weight', varargin{:});
        end
        
        function non_diag_weight = XWeight(obj, varargin)

            non_diag_weight =  obj.Categ.f_access('XWeight', varargin{:});
        end
        
        function order = f_order(obj)
            order = 0;
            error('Not known where this is used - fix function, check code');
            % order = cellfun(@(x) x(1), {obj.Category(:).Name});
        end
        
        function matrix = f_matrix(obj, vector_block_inp)
            % matrix = f_matrix(obj, vector_block_inp)
            %
            % creates Matrix A (in Energy_Structure) according to Vector Structure
            if iscell(vector_block_inp)   %input is block_type
                if sum(cellfun(@(x) numel(x),vector_block_inp)) ~= obj.Numel_Vector_Elements, error('number of elements in block does not fit to elements according so table structure'); end
                
                %NOTE: blkdiag variante - may be slow...
                %matrix = sparse(blkdiag(vector_block_inp{:}))
                
                %NOTE: via vector, function called twice
                matrix = obj.f_matrix(obj.f_vector(vector_block_inp));
                
                %NOTE: this version works via blocksize and generating the index
%                 block_size = cellfun(@(x) size(x,1),vector_block_inp);
%                 % already known from Block_Size - should be equivalent
%                 
%                 index_add = num2cell(cumsum([0;block_size(1:end-1)]));
%                 
%                 indexmesh = arrayfun(@(x) meshgrid(1:x), block_size, 'Uniformoutput', false);
%                 sparse_indices = cell2mat(cellfun(@(x,y) y + [reshape(x.',[],1), x(:)], indexmesh,index_add, 'Uniformoutput', false));
%                 
%                 matrix = sparse(sparse_indices(:,1),sparse_indices(:,2),obj.f_vector(vector_block_inp));
                
            else  %input is vector type
                if numel(vector_block_inp) ~= obj.Numel_Vector_Elements, error('Number of elements of input vector does not fit elements used in this scheme'); end
                %Note: use matrix_index function, this is fastest sparse method
                matrix_index = obj.Matrix_Index();
                matrix = sparse(matrix_index(:,1), matrix_index(:,2), vector_block_inp, obj.Numel_Datasets, obj.Numel_Datasets);
               
            end
        end
        
        
        function block = f_block(obj,vector_matrix_inp)
            % block = f_block(obj,vector_matrix_inp)
            % 
            % transforms a vector or matrix A according to Table scheme (order and filter?) to a block_structure (cell array)
            
            % if input is matrix, transform to vector first
            if issparse(vector_matrix_inp) && all(size(vector_matrix_inp) == [obj.Numel_Datasets, obj.Numel_Datasets])
                vector = sparse(f_vector(obj,vector_matrix_inp));
            else
                vector = sparse(vector_matrix_inp);
            end
            temp_cell = mat2cell(vector, obj.Block_Size.^2,1);
            
            block = cellfun(@(x,y) reshape(x, y,[]), temp_cell,num2cell(obj.Block_Size), 'uniformoutput', false);
            
            
        end
        
        function vector = f_vector(obj, matrix_block_inp)
            % vector = f_vector(obj, matrix_block_inp)
            % 
            % creates a vector from sparse A according to Table scheme
            
            if iscell(matrix_block_inp) % input is block
                vector = cell2mat(cellfun(@(x) x(:), matrix_block_inp, 'uniformoutput', false));
            else
                matrix_index = obj.Matrix_Index();
                vector = matrix_block_inp(sub2ind([obj.Numel_Datasets, obj.Numel_Datasets], matrix_index(:,1), matrix_index(:,2)));
                
            end
        end
        
        
        function [Tab] = s_print(obj, varargin)
            % s_print(obj, varargin)
            % shows Table data, Filtering possible
            
            Tab =  obj.s_table(varargin{:});
            if ~isdeployed()
                disp(['Order: ', obj.Order, ', numel: ', num2str(obj.Numel_Vector_Elements)])
            end
        end
        
        function [Tab] = s_table(obj,varargin)
             Tab = obj.Categ.s_data(varargin{:});
        end
        
        function new_tab = f_cut_filter(obj,varargin)
            %NOTE: cuts by blocks
            %NOTE: used for energy cut in born_markov function
            %TODO: check where cutting is best and what should be cut?
            %TODO: dynamical calculated objects must now be stored fixed!!!
            %(Blocksize arguments cannot be used to reconstruct for instance 
            % List Index)
            
            
            block_index = obj.Block_Index(varargin{:});  
          
            %rebuild Table:
            %obj.N_part = obj.N_part(block_index);
            %if ~isempty(obj.Spin)
            %    obj.Spin = obj.Spin(block_index);
            %end
            %obj.Energies = obj.Energies(block_index);
            %obj.Block_Size = obj.Block_Size(block_index);
            %TODO: check it - it is also in Category!!!!
            %obj.Block_Index = obj.Block_Index(block_index);
            %obj.Vector_Index = obj.Vector_Index(block_index);
            %obj.Diag_Index = obj.Diag_Index(block_index);
            
            
            %TODO: question which states shall be stored, a new, shortened Vector index
            % or the old vector indices. But since the new navigation will root on the new cut table
            % everything will be based on this new representation, an old indexing scheme might not be useful...
            % What happens with old not deriveable indices? -> Well the vector indices should always be derivable 
            % since they root on the blocks.
            
            %TODO: Why would one want to store a vector, diagonal or matrix index based on a table?
            % for calculation reasons? maybe dynamical storage? where are those indices used?
            
            % NOTE: Reduced Table is small and recalculation of Diag Indices is not so expensive.
            % It's important to delete the old ones or recompute them.
            if ~isempty(obj.Vec_Index)
                %NOTE: We may not need Vec_Index to be stored - always derived automatically
                %obj.Vec_Index = obj.Vec_Index(block_index);
            end
            obj.Vec_Index = [];
            if ~isempty(obj.Diag_Index)
                  %obj.Diag_Index = obj.Diag_Index(block_index);
            end
            obj.Diag_Index = [];
            
            
            
            
            
            %different behaviour for List_indices, since the original ones have to be restored
            
            if ~obj.Is_Cut
                additive = cumsum(obj.Block_Size) - obj.Block_Size;
                obj.Li_Index = arrayfun(@(x,y) (1:x).' + y, obj.Block_Size(block_index), additive(block_index), 'UniformOutput', false);
            else
                obj.Li_Index = obj.Li_Index(block_index);
            end
            obj.Is_Cut = true;

            
            old_indices = block_index;
            
            for n_Cat = 1:obj.Categ.Numel_Data_Types
                if obj.Categ.Data(n_Cat).Shortname == 'I'
                    obj.Categ.Data(n_Cat).Values = (1:numel(block_index)).';
                else
                    obj.Categ.Data(n_Cat).Values = obj.Categ.Data(n_Cat).Values(block_index);
                end
                obj.Categ.Data(n_Cat).Numel = numel(block_index);
                %TODO: think also about adopting variance
                %TODO: think about handling this function as function in Category!!!
                
            end
            
            obj.Categ.Numel_Data_Sets = numel(block_index);
        
            obj.Numel_Blocks = numel(block_index);
            obj.Numel_Datasets = sum(obj.Block_Size());
            obj.Numel_Vector_Elements = sum(obj.Block_Size().^2);
            new_tab = obj;
            new_tab = new_tab.f_add_categ({'Old_Indices', 'O'}, old_indices);
            

        end
        
        function new_tab = f_cut_list_index(obj,List_Index_used)
            if obj.Is_Cut
                error('only implemented for uncut tables yet!')
            end
            %NOTE: cuts by blocks, List_Index_used is a logical mask for the entries 
            %NOTE: used for energy cut in born_markov function
            %TODO: check where cutting is best and what should be cut?
            %TODO: dynamical calculated objects must now be stored fixed!!!
            %(Blocksize arguments cannot be used to reconstruct for instance 
            % List Index)
            
            block_size = cellfun(@(x) sum(x), mat2cell(List_Index_used,obj.Block_Size,1));
            block_index = find(block_size>0);            
            block_size(block_size == 0) = [];

            if ~obj.Is_Cut
                obj.Li_Index = mat2cell(find(List_Index_used), block_size,1);
            else
                %NOTE: think about it (cut tables)
            end
            obj.Is_Cut = true;
            obj.Vec_Index = [];
            
            obj.Diag_Index = [];
            
            
            %different behaviour for List_indices, since the original ones have to be restored
            
            

            
            old_indices = block_index;
            
            for n_Cat = 1:obj.Categ.Numel_Data_Types
                if obj.Categ.Data(n_Cat).Shortname == 'I'
                    obj.Categ.Data(n_Cat).Values = (1:numel(block_index)).';
                elseif obj.Categ.Data(n_Cat).Shortname == 'B'
                    obj.Categ.Data(n_Cat).Values = block_size;
                else
                    obj.Categ.Data(n_Cat).Values = obj.Categ.Data(n_Cat).Values(block_index);
                end
                obj.Categ.Data(n_Cat).Numel = numel(block_index);
                %TODO: think also about adopting variance
                %TODO: think about handling this function as function in Category!!!
                
            end
            
            obj.Categ.Numel_Data_Sets = numel(block_index);
        
            obj.Numel_Blocks = numel(block_index);
            obj.Numel_Datasets = sum(obj.Block_Size());
            obj.Numel_Vector_Elements = sum(obj.Block_Size().^2);
            new_tab = obj;
            new_tab = new_tab.f_add_categ({'Old_Indices', 'O'}, old_indices);
            

        end
        
        
        function block_index = f_vector_matrix2block_index(obj, index)
            % block_number = f_vector_matrix2block_index(obj, index)
            %
            % returns the block_number for each given vector / matrix index
            block_size = obj.Block_Size;
            if size(index,2) == 2 %matrix index
                block_end = cumsum(block_size);
                block_begin = block_end-block_size+1;
               
                row_block = arrayfun(@(x) find(x >= block_begin & x <= block_end), index(:,1));
                column_block = arrayfun(@(x) find(x >= block_begin & x <= block_end), index(:,2));
                % check, if index is in any block
                if ~all(row_block == column_block) 
                    disp('WARNING: some matrix_indices are not part of any block - they were neglected');
                end
                block_index = zeros(size(index,1),1);
                L = row_block == column_block;
                block_index(L) = row_block(L);
                
            else
                
                % make it possible for the empty Vec_Index situation:
                %arrayfun(@(x) cumsum(obj.Block_Size.^2)
                tt = cumsum(block_size.^2) - block_size.^2 ;
                if any(index < 1 | index > obj.Numel_Vector_Elements)
                    error('Error: index out of range, too low or too high')
                end
                block_index = arrayfun(@(x) find(x>tt, 1,'last'), index, 'uniformoutput', true);
                
                
            end
            
        end
        
        
        function index = f_convert_vector_matrix_index(obj,matrix_vector_index)
            
            block_size = obj.Block_Size;
            
            if size(matrix_vector_index,2) == 2 % matrix index
                % find block
                blocks = obj.f_vector_matrix2block_index(matrix_vector_index);
                
                
                % find vector_index
                L = blocks > 0;
                matrix_vector_index = matrix_vector_index(L,:);
                blocks = blocks(L);
                
                block_begin = cumsum(block_size) - block_size;
                block_internal_index = arrayfun(@(x,y,z) sub2ind([x,x], y ,z), obj.Block_Size(blocks),...
                    matrix_vector_index(:,1) - block_begin(blocks), matrix_vector_index(:,2) - block_begin(blocks));
                
                block_vector_begin = cumsum(block_size.^2) - block_size.^2;
                index = block_vector_begin(blocks) + block_internal_index;
%                index = cellfun(@(x,y) x(y), obj.Vector_Index(blocks), num2cell(block_internal_index));
            else
                blocks = obj.f_vector_matrix2block_index(matrix_vector_index);
                
                %find position in blocks and add them to 
                block_begin = cumsum(block_size) - block_size;
                block_vector_begin = cumsum(block_size.^2) - block_size.^2;
                [row, column] =  arrayfun(@(a,b) ind2sub([a,a], b) , block_size(blocks),  ...
                    matrix_vector_index - block_vector_begin(blocks));
                   
                
                index  = [row, column] + [block_begin(blocks), block_begin(blocks)];
            end
        end
        
        function start_end = f_block_start_end(obj)
            start_end = [[0; cumsum(obj.Block_Size(1:end-1))]+1, cumsum(obj.Block_Size), obj.Block_Size] ;
            warning('rethink function! where is it used?')
            
        end
            
        %DONE: clarify which way is best - so whether to use this function with the block_number or with spin /
        %particle sector; ANSW: use varargin and filter_inp to cope with both
        %TODO: define a squared index and a single index (which blocks are where) -> check if already done with
        %vector_index and index
        
        % INDEX PLAN:
        %PLAN 1: all g_..._index functions use filter_inp (varargin)
        %PLAN 2: all f_..._from_..._index translate indexing schemes
        
        function index = Vector_Index(obj,varargin)
            %warning('Here Vector_Index was used')
            %NOTE: used in born_markov_full_dynamic (209 and 523)
            block_number = obj.Block_Index(varargin{:});
             %TODO: calculate on demand, calculate all if Table is cut!? 
             % Is storage of Vec_Index necessary at all?
            if isempty(obj.Vec_Index) 
                index = [];
                for k= 1:numel(block_number) % slow implementation - Think of something faster!!!
                    index = [index; sum(obj.Block_Size(1:(block_number(k)-1)).^2) + (1:obj.Block_Size(block_number(k)).^2).'];
                end
            else %calculate dynamically
                index = cell2mat(obj.Vec_Index(block_number,1));
            end
        end
        
        function [matrix_index] = Matrix_Index(obj,varargin)

            %warning('Here Matrix_Index was used')
            %input variables are: obj.Block_Size and a presummand
            %filter according to Blocks getting the relevant block indices
            
            %TODO: think about having cut tables: where to store the Matrix_Indices?
            % calculation via the Block_Size might not work!
            index = obj.Block_Index(varargin{:});
            
            block_size = obj.Block_Size(index);
            presummand_temp = cumsum(obj.Block_Size)- obj.Block_Size;
            presummand = presummand_temp(index);
            
            
            matrix_index = cell2mat(arrayfun(@(x,y)  [reshape(meshgrid(y+(1:x)).',[],1), ...
                    reshape(meshgrid(y+(1:x)), [],1)], block_size, presummand, 'UniformOutput', false));

        end
        
        function diag_index = Diagonal_Index(obj,varargin)
            index = obj.Block_Index(varargin{:});
            %warning('Here Diagonal_Index was used')
            %NOTE: used in born_markov_full_dynamic (610, 631)
            %TODO: calculate on demand, calculate all if Table is cut!? - always dynamical derivation possible
            %(based on this (cut) Table)
            if isempty(obj.Diag_Index) 
                %TODO: calculate dynamically: find faster method?
                deg2 = cumsum([0; obj.Block_Size(1:end-1).^2]);
                diag_index_temp = arrayfun(@(x,y) (y + (1:(x+1):(x^2))'),obj.Block_Size(index), deg2(index), 'UniformOutput', false);
                diag_index = cell2mat(diag_index_temp);
            else %not needed?
                diag_index = cell2mat(obj.Diag_Index(index));
            end
        end
        
        function Row_Index(obj,varargin)
            error('use List_Index instead')
        end
        
        function list_index = List_Index(obj,varargin)
           
            %TODO: think about having cut tables: where to store the List_indices?
            %NOTE: solution: function Old_List_Index -> Li_Index; List_Index -> fresh calculated
            % calculation via the Block_Size might not work!
            index = obj.Block_Index(varargin{:});
            if ~isempty(obj.Li_Index) && ~obj.Is_Cut
                %TODO.... Think about what happens with cut table
                list_index = cell2mat(obj.Li_Index(index));
            else
            
                if nargin < 2 || isempty(varargin)
                    list_index = 1:obj.Numel_Datasets;
                else
                    
                    additive = cumsum(obj.Block_Size) - obj.Block_Size;
                    list_index = cell2mat(arrayfun(@(x,y) (1:x).' + y, obj.Block_Size(index), additive(index), 'UniformOutput', false));
                end
            end
        end
        
        function list_index = Old_List_Index(obj,varargin)
            
            
            index = obj.Block_Index(varargin{:});
            if ~isempty(obj.Li_Index) && obj.Is_Cut
                %TODO.... Think about what happens with cut table
                list_index = cell2mat(obj.Li_Index(index));
            else
                
                if nargin < 2 || isempty(varargin)
                    list_index = 1:obj.Numel_Datasets;
                else
                    
                    additive = cumsum(obj.Block_Size) - obj.Block_Size;
                    list_index = cell2mat(arrayfun(@(x,y) (1:x).' + y, obj.Block_Size(index), additive(index), 'UniformOutput', false));
                end
            end
        end
        
        
        function K2 = f_choi(obj, K)
            
            [i1,j1,s1] = find(K);
            block_size = obj.Block_Size;
            
            % instead of
            %
            %A = [obj.f_convert_vector_matrix_index(i), obj.f_convert_vector_matrix_index(j)];
            % using:
            %blocks = obj.f_vector_matrix2block_index(i1);
            
            block_vector_begin = cumsum(block_size.^2) - block_size.^2 ;
            block_begin = cumsum(block_size) - block_size;
            
            %superrow i
            blocks_i = arrayfun(@(x) find(x>block_vector_begin, 1,'last'), i1, 'uniformoutput', true);
                      
            [row_i, column_i] =  arrayfun(@(a,b) ind2sub([a,a], b) , block_size(blocks_i), i1 - block_vector_begin(blocks_i));
            
            %supercolumn j
            blocks_j = arrayfun(@(x) find(x>block_vector_begin, 1,'last'), j1, 'uniformoutput', true);
                
            [row_j, column_j] =  arrayfun(@(a,b) ind2sub([a,a], b) , block_size(blocks_j), j1 - block_vector_begin(blocks_j));
             
            % superindices (a,b,c,d)
            A = [row_i + block_begin(blocks_i), column_i + block_begin(blocks_i), row_j + block_begin(blocks_j), column_j + block_begin(blocks_j)];
            
            
            in1 = sub2ind([obj.Numel_Datasets, obj.Numel_Datasets], A(:,1), A(:,3));
            in2 = sub2ind([obj.Numel_Datasets, obj.Numel_Datasets], A(:,2), A(:,4));
            
            %cutting zero rows and columns
            [~,~,rank1] = unique(in1);
            [~,~,rank2] = unique(in2);
            K2 = sparse(rank1,rank2,s1); 
            
            
            
           
   
        end
        
        function [choi_herm, choi_trace, choi_pos_ev, choi_ev] = f_choi_analyse(obj, K)
            % evaluates hermiticity of Choi matrix, trace and eigenvalues
            K2 = obj.f_choi(K);
            
            % =============  hermiticity: =================
            T = 1/2*(K2 - K2');
            %NOTE: Divide by nnz? -> NO! Definition as in ERMEA Paper
            choi_herm = sqrt(sum(abs(T(:).^2)));
            
            % ================  trace: ====================
            %NOTE: trace preserving condition of dynamical map not checked -> do outside
            choi_trace = trace(K2);
            
            if nargout > 2
                % =============  eigenvalues: =================
                % block matrix sorting
                p = symamd(K2);
                K3 = K2(p,p);
                
                %recognizing blocks
                [i2,j2,s2] = find(K3);
                %spy(K3)
                
                
                %NOTE: presort with symamd - then according to the sparse ordering a block on a matrix with a non-zero
                %diagonal is found if two diagonal elements occur subsequently, thus the border of a block is found,
                %cumsum generates the block number
                diagonal_match = i2 == j2;
                block_id = cumsum([1;diagonal_match(1:end-1)] == diagonal_match & diagonal_match == 1);
                numel_block_elements = hist(block_id, 1:max(block_id)).';
                blocks = arrayfun(@(x,y) sparse(i2(x+(1:y))-i2(x+1)+1, j2(x+(1:y)) - j2(x+1) + 1, s2(x+(1:y))), cumsum([0;numel_block_elements(1:end-1)]), numel_block_elements, 'uniformoutput', false);
                %NOTE: control blocks
                %AA = blkdiag(blocks{:});
                
                choi_ev = cell2mat(cellfun(@(x) eig(full(x)), blocks, 'uniformoutput', false));
                
                choi_pos_ev = sum(abs(choi_ev) - choi_ev);
            end
            
%             if length(K) < 1500
%                 [xx,yy] = meshgrid(1:obj.Numel_Vector_Elements);
%                 %[xx,yy] = meshgrid(188:195);
%                 A = [obj.f_convert_vector_matrix_index(yy(:)), obj.f_convert_vector_matrix_index(xx(:))];
%                 % check for Hermiticity
%                 
%                 in_1 = sub2ind([obj.Numel_Datasets, obj.Numel_Datasets], A(:,2), A(:,1));
%                 in_2 = sub2ind([obj.Nume_Datasets, obj.Numel_Datasets], A(:,4), A(:,3));                
%                 KK2 = sparse(in_1, in_2, K(:), obj.Numel_Datasets.^2, obj.Numel_Datasets.^2);
%                 rred1 = sum(KK2~=0,1);
%                 rred2 = sum(KK2~=0,2);
%                 LL = (rred1==0) & (rred2'==0);
%                 KK3 = KK2(~LL,~LL);
%                 
%                 
%                 
%                 in1 = sub2ind([obj.Numel_Datasets, obj.Numel_Datasets], A(:,1), A(:,3));
%                 in2 = sub2ind([obj.Numel_Datasets, obj.Numel_Datasets], A(:,2), A(:,4));
%                 
%                 K2 = sparse(in1, in2, K(:), obj.Numel_Datasets.^2, obj.Numel_Datasets.^2);
%                 
%                 red1 = sum(K2~=0,1);
%                 red2 = sum(K2~=0,2);
%                 L = (red1==0) & (red2'==0);
%                 K3 = K2(~L,~L);
%                 
%                 
%                 T = K3 - K3';
%                 %NOTE: Divide by nnz?
%                 choi_herm = sqrt(sum(abs(T(:).^2)));
%                 
%                 K3_diag = K3 - diag(diag(K3));
%                 red_1 = sum(K3_diag~=0,1);
%                 red_2 = sum(K3_diag~=0,2);
%                 L2 = (red_1==0) & (red_2' == 0);
%                 
%                 K4 = K3(~L2,~L2);
%                 ew1 = diag(K3);
%                 ew1 = ew1(L2);
%                 %develop a rough eigenvalue estimation!
%                 
%                 if length(K4) <1000
%                     choi_ew = eig(full(K4));
%                 else
%                     choi_ew = nan(length(K4),1);
%                 end
%                 
%                 choi_ew = [ew1; choi_ew];
%                 
%             else
%                 warning('K larger than 1000 - no choi number calculated!')
%                 K2 = [];
%                 choi_herm = inf;
%                 choi_ew = inf;
%             end
        end
        
        
        function [super_index, K2] = f_choi_index_permutation(obj, K)
            warning('OLD - check where used!')
            [xx,yy] = meshgrid(1:obj.Numel_Vector_Elements);
            %[xx,yy] = meshgrid(188:195);
            A = [obj.f_convert_vector_matrix_index(xx(:)), obj.f_convert_vector_matrix_index(yy(:))];
            
            A1 = A(:, [1,3])
            in1 = sub2ind([obj.Numel_Datasets, obj.Numel_Datasets], A(:,1), A(:,3));
            in2 = sub2ind([obj.Numel_Datasets, obj.Numel_Datasets], A(:,2), A(:,4));
            
            K2 = sparse(in1, in2, K(:), obj.Numel_Datasets.^2, obj.Numel_Datasets.^2);
            
            
            
            %should return the same!!!
            Vector_Index_1 = obj.f_convert_vector_matrix_index(A(:,[1,2]));
            Vector_Index_2 = obj.f_convert_vector_matrix_index(A(:,[3,4]));
            super_index = sub2ind([obj.Numel_Vector_Elements, obj.Numel_Vector_Elements], Vector_Index_1, Vector_Index_2);
            
        end
        
        function sigma_vector = f_create_sigma_vector(obj, type)
            switch lower(type(1))
                case 'r' % random
                    block_weight = rand(obj.Numel_Blocks,1);
                    block_weight = num2cell(block_weight/sum(block_weight));
                    temp = arrayfun(@(x) rand(x)-0.5 + 1i * (rand(x)-0.5),obj.Block_Size, 'Uniformoutput', false);
                    temp = cellfun(@(x,y) (x*x')/sum(abs(x(:)).^2) * y, temp, block_weight, 'uniformoutput', false);
                     
                    sigma_vector = obj.f_vector(temp);

                    
                case 'o' % ones = maximal entangled state = left vacuum
                    sigma_vector = obj.f_vector(eye(obj.Numel_Datasets)/obj.Numel_Datasets);
                case 'g' % ground state
                    temp = sparse(obj.Numel_Datasets, obj.Numel_Datasets);
                    if ismember('E', obj.Order)
                        min_en = min(obj.Energies);
                        row_index = obj.List_Index('E', min_en);
                        
                    elseif  ismember('A', obj.Order)
                        min_en = min(obj.Average_Energies);
                        row_index = obj.List_Index('A', min_en);
                    end
                    temp(row_index, row_index) = eye(numel(row_index))/numel(row_index);
                    
                    sigma_vector = obj.f_vector(temp);
            
                    %temp(obj.List_Index('A', 
            end
            
        end
        
        
         %NOTE: To be deleted
        function block_index = g_block_index(obj, varargin)
            error('use Index instead');
           
        end
        
      
        % TODO: ------------------ check, if needed
        function new_index = f_index_translation(obj, index_in, table)
            %get a new index according to different Energy / Block scheme (Energies sorted uniquely!!!
            warning('Warning: Old????')
            vector = sparse(1,obj.Numel_Vector_Elements);
            vector(index_in) = 1:numel(index_in); % to keep order
            full = obj.f_matrix(vector);
            new_vector = table.f_vector(full);
            raw_index = find(new_vector);
            [~, new_order] = sort(new_vector(new_vector > 0));
            new_index = raw_index(new_order);
            

        end
        
        function Sig = f_eigenvector_to_rho(obj,rho)
            sigma_tolerance = 10^-10;
            Sig_vec_denoise = denoise(rho, sigma_tolerance);
            
            Values = obj.f_block(Sig_vec_denoise);
            
            block_weight = cellfun(@(x) sum(abs(full(x(:)))), Values);
            L_ind = find(block_weight > sigma_tolerance);
            
            weight = cell2mat(cellfun(@(x,y) repmat(sum(abs(full(x(:)))),y,1), Values(L_ind), num2cell(obj.Block_Size(L_ind)), 'uniformoutput', false));
            argument = full(cell2mat(cellfun(@(x) atan2(imag(diag(x)), real(diag(x))), Values(L_ind), 'uniformoutput', false)));
            argument = mod(argument, 2*pi);
            if mean(argument) < pi % in order to avoid jumps!!!
                argument(argument > 5.5) = argument(argument > 5.5) - 2*pi;
            else
                argument(argument < 1) = argument(argument < 1) + 2*pi;
            end
            
            % if matrix is truely hermitian then only two values are possible!
            [unique_phase,~,frequency_of_phase] = uniquetol(argument(weight> sigma_tolerance), 10^-2);
            phase_correction = median(unique_phase(frequency_of_phase));
            
            if numel(unique_phase) > 2%var(argument, weight)^2 > sigma_tolerance
                warning('density matrix has not the same complex phase on all diagonal elements!')
                %NOTE: density matrix should be Hermitian, but could have negative
                %eigenvalues - therefore discrimination
            end
            
            if numel(unique_phase) == 2
                warning(['Density matrix is negative!!!'])
            end
            
            Sig_vec_denoise = denoise(Sig_vec_denoise * exp(-1i*phase_correction), sigma_tolerance);
            
            S = obj.f_matrix(Sig_vec_denoise);
            % NOTE: trace norm always 1 - not equal to trace(abs(rho))!!!
            trace_abs = trace((S'*S)^(1/2));
            Sig_vec_denoise = denoise(Sig_vec_denoise./trace_abs, sigma_tolerance);
            Sig = Sigma([], Sig_vec_denoise, obj);
        end
        
        
        function [operator_norm, L_2_norm, L_spec_norm] = f_operator_norm(obj, K, test_rho)
            if nargin == 3
                flag_special_rho = true;
                n_rho = size(test_rho,2);
            else
                flag_special_rho = false;
                n_rho = 4000;
            end
         
%             for ind = 1:obj.Numel_Vector_Elements
%                 T = zeros(obj.Numel_Vector_Elements,1);
%                 T(ind) = 1; %obj.Block_Size(obj.f_vector_matrix2block_index(index));
%                 rho{ind} = obj.f_matrix(T) + obj.f_matrix(T)' + eye(obj.Numel_Datasets); %/obj.Block_Size(obj.f_vector_matrix2block_index(index))
%                 %eig(rho{ind})
%                
%             end
            
            


            for k = 1:n_rho %obj.Numel_Vector_Elements
                
                if ~flag_special_rho
                    
                    % --- special chosen basis
                    %X = Sigma([],sparse(rho{k}),obj);
                    
                    % --- random Hermitian matrices
                    %X = rand(obj.Numel_Datasets) - 0.5 + 1i*(rand(obj.Numel_Datasets)-0.5);
                    %X = Sigma([],sparse(X + X'),obj);
                    
                    % --- random cut matrices
                    %X = Sigma([],sparse(random_rho(obj.Numel_Datasets)),obj);
                    
                    
                    % --- random density matrices
                    X = Sigma([],[],obj);
                else
                    X = obj.f_eigenvector_to_rho(test_rho(:,k))
                    
                end
                A = X.f_vector;
                %A = obj.f_vector(random_rho(obj.Numel_Datasets));
                B = obj.f_matrix(A);
                T = obj.f_matrix(K*A);

                operator_norm(k) = trace( (T' *T)^(1/2))./trace( (B'*B)^(1/2));
                
                %trace(abs(obj.f_matrix(K * A)))./trace(abs(B));
                L_2_norm(k) = sqrt(sum(sum(abs(K * A).^2)))./sqrt(sum(abs(B(:)).^2));
                L_spec_norm(k) = max(eig(full(obj.f_matrix(K*A))))./(max(eig(full(B))));
                %if operator_norm(k) > 1.001
                %   disp(X.s_dom_blocks)
                %    disp(operator_norm(k))
                %end
            end
            
            
            
            
            a = hist(real(operator_norm), 0.98:0.01:1.1)
            L_2_norm = max(L_2_norm);
            L_spec_norm = max(L_spec_norm);
                        
            operator_norm = max(operator_norm);
%             
%             t = linspace(-1,1,100);
%             for k=1:numel(t)
%                 for j = 1:numel(t)
%                     R = zeros(3);
%                     R(1,1) = 0.4;
%                     R(2,2) = 0.1;
%                     R(3,3) = 0.2;
%                     R(1,3) = t(k);
%                     R(2,3) = t(j);
%                     T = R + R';
%                     ew = eig(T);
%                     E(k,j) = sum(ew(ew<0));
%                 end
%             end
%             E(E>-10^-5) = 1;
%             [tx,ty] = meshgrid(t);
%             pcolor(tx,ty,E)
            
            
            
          
            
        end
                
                
        
        
        function [fig1] = s_plot(obj,varargin)
            
            block_index = obj.Block_Index(varargin{:});
            a = 0.5;
            hold on
            
            xtick = zeros(1,numel(block_index));
            for n_block = 1:numel(block_index)
                
                blocksize = obj.Block_Size(block_index(n_block));
                
                plot(a + blocksize*[0,1,1,0,0],a + blocksize * [0,0,1,1,0], 'g-', 'linewidth', 1.5)
                
                xtick(n_block) = a + blocksize/2;
                
                a = a + blocksize;
                
            end
            %TODO: workaround if energy or weight is not known
            
            
            if ismember('A', obj.Order)
                energies = obj.Average_Energies(block_index);
            elseif ismember('E', obj.Order)
                energies = obj.Energies(block_index);
            else
                energies = '';
            end
            
            
            part = obj.N_part(block_index);
            if ismember('S', obj.Order)
                spin = obj.Spin(block_index);
                part_spin = arrayfun(@(x,y) ['N',  num2str(x), ', S', num2str(y)], part, spin, 'uniformoutput', false);
            else
                part_spin = arrayfun(@(x) ['N',  num2str(x)], part, 'uniformoutput', false);
            end
            
            
            char_array_en = num2str(energies,4);
            if ismember('W', obj.Categ.f_shortnames)
                char_array_weight = num2str(obj.Weight(block_index),4);
            else
                char_array_weight = num2str('');
            end
            labels_weight = mat2cell(char_array_weight, ones(length(block_index),1), size(char_array_weight,2));
            set(gca, 'xtick', xtick)
            set(gca, 'xticklabel', part_spin)
            set(gca, 'xTickLabelRotation', 90)
            xlabel('particle and spin number')
            if ~isempty(energies)
                
                labels_energies = mat2cell(char_array_en, ones(length(block_index),1), size(char_array_en,2));
                set(gca, 'yTick',xtick)
                set(gca,'yTickLabel', labels_energies)
                ylabel('energies')
                set(gca, 'YDir', 'reverse')
                
                
            end
            yyaxis right
            
             a = 0.5;
            hold on
            
            xtick = zeros(1,numel(block_index));
            for n_block = 1:numel(block_index)
                
                blocksize = obj.Block_Size(block_index(n_block));
                
                plot(a + blocksize*[0,1,1,0,0],a + blocksize * [0,0,1,1,0], 'g-', 'linewidth', 1.5)
                
                xtick(n_block) = a + blocksize/2;
                
                a = a + blocksize;
                
            end
            set(gca, 'YDir', 'reverse')
            set(gca, 'yTick',xtick)
            set(gca,'yTickLabel', labels_weight)
            ylim([0.5,a])
            yyaxis left
            ylim([0.5,a])


            grid on

            
        end
        
        function [fig1] = s_plot_colorwidth(obj,varargin)
            
            block_index = obj.Block_Index(varargin{3:end});
            a = 0.5;
            hold on
            
            xtick = zeros(1,numel(block_index));
            for n_block = 1:numel(block_index)
                
                blocksize = obj.Block_Size(block_index(n_block));
                
                plot(a + blocksize*[0,1,1,0,0],a + blocksize * [0,0,1,1,0], 'color', varargin{1}, 'linewidth', varargin{2}, 'linestyle', '-', 'marker', 'none')
                
                xtick(n_block) = a + blocksize/2;
                
                a = a + blocksize;
                
            end
            %TODO: workaround if energy or weight is not known
            
            
            if ismember('A', obj.Order)
                energies = obj.Average_Energies(block_index);
            elseif ismember('E', obj.Order)
                energies = obj.Energies(block_index);
            else
                energies = '';
            end
            
            
            part = obj.N_part(block_index);
            if ismember('S', obj.Order)
                spin = obj.Spin(block_index);
                part_spin = arrayfun(@(x,y) ['N',  num2str(x), ', S', num2str(y)], part, spin, 'uniformoutput', false);
            else
                part_spin = arrayfun(@(x) ['N',  num2str(x)], part, 'uniformoutput', false);
            end
            
            
            char_array_en = num2str(energies,4);
            if ismember('W', obj.Categ.f_shortnames)
                char_array_weight = num2str(obj.Weight(block_index),4);
            else
                char_array_weight = repmat(num2str(' '),length(block_index),1);
            end
            labels_weight = mat2cell(char_array_weight, ones(length(block_index),1), size(char_array_weight,2));
            set(gca, 'xtick', xtick)
            set(gca, 'xticklabel', part_spin)
            set(gca, 'xTickLabelRotation', 90)
            xlabel('particle and spin number')
            if ~isempty(energies)
                
                labels_energies = mat2cell(char_array_en, ones(length(block_index),1), size(char_array_en,2));
                set(gca, 'yTick',xtick)
                set(gca,'yTickLabel', labels_energies)
                ylabel('Average energies', 'color','r')
                set(gca, 'YDir', 'reverse')
                ticklabels = get(gca,'YTickLabel');
                % prepend a color for each tick label
                ticklabels_new = cell(size(ticklabels));
                for i = 1:length(ticklabels)
                    ticklabels_new{i} = ['\color{red} ' ticklabels{i}];
                end
                % set the tick labels
                set(gca, 'YTickLabel', ticklabels_new);
                
            end
            yyaxis right
            
             a = 0.5;
            hold on
            
            xtick = zeros(1,numel(block_index));
            for n_block = 1:numel(block_index)
                
                blocksize = obj.Block_Size(block_index(n_block));
                
                plot(a + blocksize*[0,1,1,0,0],a + blocksize * [0,0,1,1,0], 'color', varargin{1}, 'linewidth', varargin{2}, 'linestyle', '-', 'marker', 'none')
                
                xtick(n_block) = a + blocksize/2;
                
                a = a + blocksize;
                
            end
            set(gca, 'YDir', 'reverse')
            set(gca, 'yTick',xtick)
            set(gca,'yTickLabel', labels_weight)
            ylabel('Weights')
            ylim([0.5,a])
            yyaxis left
            ylim([0.5,a])


            grid on

            
        end
        
        
        function [fig1] = s_plot_old(obj, A)
            % define a routine to visualize block structure, and to fill in
            % a) as matrix - sparse matrix
            % b) as block  - already block structure
            % c) as single value - indicating one block
            % d) as vector - vector index
            %TODO: define colors, label structure, ticks on the side?
            %TODO: rethink input - maybe an arbitrary block structure could be plotted, 
            % or only the relevant ones (rounding on matrix A?)
            max_index = obj.Numel_Vector_Elements;
            if nargin == 2
                a = 0.5;
                if issparse(A) % matrix structure (at least it should be sparse)
                    % A = A
                else
                    if ~iscell(A) && numel(A) == 1 % A is a number indicating the single block
                        indices = obj.Vector_Index(A);
                        vector = sparse(1,max_index);
                        vector(indices) = 1:numel(indices);
                        A = obj.f_matrix(vector);
                    else %block or vector structure
                        
                        
                        A = obj.f_matrix(A);
                    end
                    
                end
                
                
            else
                A = sparse(obj.Numel_Datasets, obj.Numel_Datasets);
                a = 1; %TODO: find a way to start with spy at (1,1)
            end
            
            spy(A);
            %continue by drawing the corners of the blocks
            hold on
            
            for n_block = 1:numel(obj.Block_Index)
                
                blocksize = obj.Block_Size(n_block);
                
                plot(a + blocksize*[0,1,1,0,0],a + blocksize * [0,0,1,1,0], 'g-', 'linewidth', 2)
                a = a + blocksize;
            end
            
            
            
            fig1 = gcf;
            
            
        end
        
        function [] = s_Blocks(obj)
            figure
            
            n_part = unique(obj.N_part);
            subplot_rows = min(numel(n_part),3);
            subplot_cols = ceil(numel(n_part)/subplot_rows);
            figure
            plot_settings(false)
                
            for k = 1:numel(n_part)
                n_p = n_part(k);
                block_size = obj.Block_Size('N', n_p);
                if ismember('E', obj.Order)
                    energ = obj.Energies('N', n_p);
                elseif ismember('A', obj.Order)
                    energ = obj.Average_Energies('N',n_p);
                else
                    energ = ones(size(block_size));
                end
                
                
                
                    
                subplot(subplot_rows, subplot_cols, k)
                plot(energ, block_size, 'x', 'linestyle', 'none')
                xlabel('energy')
                ylabel('block_size')
                    
                    
                
                
            end
            
        end
        
    end
    
    
end
