classdef Category 
    %author: Gerhard Dorn
    %date: august 16
    %use: store and access data (1-1 mappings)
    %     search (filter)
    %     
    %features: - 
    %          - 
    %          - 
    %functionality: - append(name, Values, Tol, Variance)
    %               - access(name, filter)
    %               - filter
    %               - shortlist
    
    properties
        Data %stores the data columns
        List %TODO: check if unique indexing is really needed, there may be problems with indexing for future extensions??
        Numel_Data_Types
        Numel_Data_Sets
        
        
    end
    
    
    properties (Access = private)
        
    end
    
    methods
        function obj = Category(Name, Values, Tol, Variance)
           obj.Data = []; % TODO: Maybe initialize with needed structure - empty value seems also to work
           % Name shall can consist of Cell: {'Name', 'Shortname'}
           obj.List = containers.Map('KeyType','char','ValueType','int32');
           obj.Numel_Data_Types = 0;
           obj.Numel_Data_Sets = 0;
           
           if nargin >= 2 %input present and valid, TODO: make an input check
              % if nargin < 5 || isempty(Type), Type = 'double'; end
               if nargin < 4 || isempty(Variance), Variance = []; end
               if nargin < 3 || isempty(Tol), Tol = -1; end
               
               obj = obj.f_append(Name, Values, Tol, Variance);
           end
        end
        
        
        function obj = f_append(obj, Name_inp, Values, Tol, Variance)
            % appends the Data by a new Category (a new data set with new name and shortname)
            % Name_inp shall can consist of Cell: {'Name', 'Shortname'}
            %if nargin < 6 || isempty(Type), Type = 'double'; end
            if nargin < 5 || isempty(Variance), Variance = []; end
            if nargin < 4 || isempty(Tol), Tol = -1; end
            % more attributes can be extended - 
            % default initialization shall be used to destinguish if existent or not
            % Name_inp is cell for dedicated shortname or not if deduced from first letter
            
            %check input
            %TODO: to be extendend
            if (~iscell(Name_inp)) || (numel(Name_inp) ~= 2)
                warning('No shortname defined explicitly, first letter will be taken')
                Shortname = Name_inp(1);
                Name = Name_inp;
            else
                if (~ischar(Name_inp{2})) || (numel(Name_inp{2}) ~= 1) 
                    error('Shortname has to consist of one letter - case sensitive and unique')
                end
                if obj.List.isKey(Name_inp{2})
                    error(['Shortname already defined: ',obj.f_shortnames])
                end
                Shortname = Name_inp{2};
                Name = Name_inp{1};
            end
            
            % new_index
            new_index = obj.Numel_Data_Types + 1;
            
            % save values
            obj.Data(new_index).Name = Name;
            
            
            
            
            obj.Data(new_index).Values = Values;
          
            obj.Data(new_index).Tol = Tol;
            obj.Data(new_index).Variance = Variance;
            obj.Data(new_index).Shortname = Shortname;
            obj.Data(new_index).Numel = numel(Values);
            %TODO: check if .Numel is used anyway!!!
            
            % update number of data sets in the initial case
            if obj.Numel_Data_Types == 0
                obj.Numel_Data_Sets = numel(Values);
            end
            
            % update number of Data types
            obj.Numel_Data_Types = numel(obj.Data);
            
            % update container maps (lists)
            obj.List(Name) = new_index;
            obj.List(Shortname) = new_index;

        end
        
        function obj = f_remove(obj, Name_inp)
            % removes a Category from the data set
            if iscellstr(Name_inp)
                Name = Name_inp{1};
            else
                Name = Name_inp;
            end
            
            
            if ~ischar(Name), error('f_remove requires the name of the category to be removed as string!'); end
            
            index = obj.List(Name);
            obj.Data(index) = [];
            obj = obj.f_recreate_lists;
            
            if iscellstr(Name_inp) && numel(Name_inp) > 1
                obj = obj.f_remove(Name_inp(2:end));
            end
             
            
        end
        
        
        function values = f_access(obj, Name, varargin)
            
            if nargin < 3 || isempty(varargin) 
                % if no filter is to be applied, return all data (values)
                values = obj.Data(obj.List(Name)).Values;
            else
                filter = obj.filter_inp(varargin{:});
                values = obj.Data(obj.List(Name)).Values(filter);
                
                %DONE: define the task of filter_inp in detail!!!
            
            
                
            end
        end
        
        function LL = filter_inp(obj, varargin)
            % redefine filter_inp function!
            % criterias:
            %   - shall be fast for single numerical values (prevent ismember function)
            %   - shall be dynamic for extended filter searches
            %   - shall be fast in terms of indexing - logical values vs. indexes, define criteria? 
            %   - shall work as if a normal variable is querried (index, subindex, logical variable)
            % the principle search algorithms return logical variables - do sparse version exist?
            % NOTE: work now with logical variables - maybe there are logical sparse vectors
       
            input = varargin(:);
            num_varargin = numel(input);
            filter_sequence = input{1};
            
            
            if ~ischar(input{1})  %NOTE: normal indexing mode, three possibilities
                % - logical vector
                % - index
                % - colon (:)
                
                if num_varargin == 1 && islogical(input{1})
                    LL = input{1};
                end
                
                %NOTE: preserves sequence!!!
                if num_varargin == 1 && isnumeric(input{1})
                    LL = input{1};
                end
                
                %NOTE:  reactivate if Category contains also square matrices -  this may be slow 
                %if num_varargin == 2  % more difficult to deal with matrices
                 
                %    LL = false(obj.Data(1).Numel, obj.Data(1).Numel);
                %    LL(input{1}, input{2}) = true;
                %end                
                
            else %NOTE: sequence lost!!!
                    
                
            
            
                
                % NOTE: error catched later when filter_sequence not found in Map
                %FILTER_CAT_LABEL = obj.f_shortnames();
                %if ~all(ismember(filter_sequence, FILTER_CAT_LABEL)),
                %    error(['Filter categories ', filter_sequence, ' do not match ', FILTER_CAT_LABEL]);
                %end
                
                num_filter = numel(filter_sequence);
                if num_filter ~= num_varargin - 1 % given filter categories 'NE' corresponds to number of filters input - the one for selection mode
                    error('number of filter categories does not match delivered filters');
                end
                
                %filter sequence
                % fixed categories, names should be unique, in filter sequence and in category labelling:
                % E ... Energy, include tolerance
                % S ... Spin, include variance
                % N ... Particle sector, include variance
                % D ... Double occupation / Degeneracy
                
                
                
                
                %used to identify if Tolerance is applicable (-1 Tolerance means Tolerance feature not used, 0
                %Tolerance can also be defined but has no influence.
                FILTER_CAT_TOL = [obj.Data(:).Tol] > 0;
                
                
                L = cell(1,num_filter);
                LL = true(obj.Numel_Data_Sets, 1);
                
                
                for k = 1:num_filter
                    if obj.List.isKey(filter_sequence(k))
                        data_index = obj.List(filter_sequence(k)); %use Map to get right Data index
                        %FILTER_CAT = obj.Data(data_index).Value;
                        
                        filter_temp = input{k+1}; %get filters
                        if isscalar(filter_temp)     %single value - fast version
                            if FILTER_CAT_TOL(data_index)
                                L{k} = (obj.Data(data_index).Values < filter_temp +  obj.Data(data_index).Tol) & (obj.Data(data_index).Values > filter_temp - obj.Data(data_index).Tol);
                            else
                                L{k} = obj.Data(data_index).Values == filter_temp;
                            end
                            
                        elseif isnumeric(filter_temp) && isvector(filter_temp)  %vector
                            if FILTER_CAT_TOL(data_index)
                                %use tolerance (only works for newer Matlab > 2015a)
                                vers = version();
                                if str2double(vers(end-3:end-2)) >= 15
                                    L{k} = ismembertol(obj.Data(data_index).Values, filter_temp, obj.Data(data_index).Tol);
                                else
                                    
                                    filter_temp = round(filter_temp, -floor(log10(obj.Data(data_ind).Tol )));
                                    L{k} = ismember(obj.Data(data_index).Values, filter_temp);
                                end
                            else
                               L{k} = ismember(obj.Data(data_index).Values, filter_temp);
                            end
                        elseif iscell(filter_temp)  %range
                            if numel(filter_temp) ~= 2, error('for range filter use cell with two entries'); end
                            L{k} = obj.Data(data_index).Values >= filter_temp{1} & obj.Data(data_index).Values <= filter_temp{2};
                            
                            
                        elseif isempty(filter_temp) %no value
                            L{k} = true;   %TODO: should no arrise
                        elseif islogical(filter_temp) && all(size(filter_temp) == size(LL)) %if input is already a logical vector
                            L{k} = filter_temp; %TODO: should not arrise
                            
                        else
                            
                        end
                        %NOTE: potential to add an OR here - new function
                        LL = LL & L{k};
                    else
                        if ~isempty(input{k+1}) && filter_sequence(k)~= '!'
                            warning(['Category|filter_inp: the Category ', filter_sequence(k), ' is not part of the categories (table_tag): ', obj.f_shortnames, ' and will be ignored!.'])
                        end
                    end
                end
            end
            
        end
                
        function shortnames = f_shortnames(obj, indices)
            %TODO: check if firstones was needed somewhere...
            shortnames = [obj.Data(1:end).Shortname];
            if nargin == 2 
                shortnames = shortnames(indices);
            end
        end
        
        function longnames = f_longnames(obj, indices)
            longnames = {obj.Data(1:end).Name};
            if nargin == 2 
                longnames = longnames(indices);
            end
        end
        
        % may not be used so often
        function obj = f_recreate_lists(obj)
            name_keys = {obj.Data(:).Name};
            short_keys = {obj.Data(:).Shortname};
            obj.Numel_Data_Types = numel(obj.Data);
            values = 1:obj.Numel_Data_Types;
            
            % reference: obj.List(Name) = new_index;
            % alternative way is to create two maps and to add them by concatenating
            obj.List = containers.Map([name_keys, short_keys], [values, values]);

            
        end
        
        
        function obj = f_add(obj, Shortnames, varargin)
            % adds a single data set, all existing categories have to be set, otherwise warning and empty (or
            % zero?) 
            % TODO: define cases where you don't fill all entries - how to deal with it...
            
            
            %check Shortname range
            %if ~isstring(Shortnames)
            if ~all(ismember(Shortnames, obj.f_shortnames)) && ~all(obj.f_shortnames, ['I',Shortnames])
                error('The dataset added to Data Container, has to have the same Categories except index')
                
            end
            %find new index
            
            data_set_index = obj.Numel_Data_Sets + 1;
            
            
            % add Index first
            obj.Data(obj.List('I')).Values(data_set_index) = data_set_index;
            obj.Data(obj.List('I')).Numel = obj.Data(obj.List('I')).Numel + 1;
            
            for k = 1:numel(Shortnames)
                type_index = obj.List(Shortnames(k));
                obj.Data(type_index).Values(data_set_index) = varargin{k};
                obj.Data(type_index).Numel = obj.Data(type_index).Numel + 1;
                
            end
            
            obj.Numel_Data_Sets = obj.Numel_Data_Sets + 1;
            
            %TODO: check speed - will not be able to adopt as repetitive procedure - think about different scheme!
            
            %TODO: maybe adopt a function for update of individual class.
            %construction of having always the same number of entries in each dataset is not ultimativly clear so
            %far - there may be other applications            
            
            
            
            %loop over all varargin sets - should be a fast procedure!
            
            
            
            
            
            
            
            
            
        end
        
        %-------------------display function-------------------------
        function  tab = s_names(obj)
            indices = (1:obj.Numel_Data_Types)';
            names = {obj.Data(:).Name};
            shortnames = {obj.Data(:).Shortname};
            for k = 1:obj.Numel_Data_Types
                temp = obj.Data(k).Values;
                s(k,1) = whos('temp');
            end
            
            
            tab = table(indices(:), names(:), shortnames(:), {s.class}', [s.bytes]', 'variablenames', {'Index', 'Name', 'Shortname', 'Type', 'Bytes'});
            disp(tab)
            
        end
            
        function tab = s_data(obj, varargin)
           
                      
            names = obj.f_longnames;
          
            tab = array2table([obj.Data(:).Values], 'VariableNames', names); 
            
            
            if nargin > 1 && ~isempty(varargin)
                LL = filter_inp(obj, varargin{:});
                tab = tab(LL,:);
            end
            
            %disp(obj.s_filter(varargin{:}));
           
        end
        
        
        function string = s_filter(obj, varargin)
            %TODO: show filter
            
            
        end
        
    end
    
    
end