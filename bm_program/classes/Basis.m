classdef Basis < handle
    %author: Gerhard Dorn
    %date: august 15
    %use: create a fermionic basis and save it as an integer number
    %features: - define number of sites (N_site), number of electrons (Np_input), total spin (Sp_input) (1 particle = spin of +-1)
    %          - restrict double / zero occupation (Docc)
    %          - spinless fermions (Num_of_Spin)
    %functionality: - easy access spin sectors, address all properties
    
    %TO DO: make some properties private so that tab shows just relevant
    %commands / functions
    
    properties
        N_site          % Number of system sites
        N_spin          % - Number of spins used (standard is 2)     
        Smaller_Binary_Sector_Used_For % indicates which spin is on the right side of binary 01011|00100 = 356
        Numel_Bas       % number of basis elements
        %Bas_Dez         % list of all basis states converted into dezimal representation (sorted, by part. num. | spin | increasing number)
        %N_part          % saves a list of the particle number of all basis states
        %Spin            % saves a list of the netto spin of all basis states (spin 1 particles)
        %Double_Occ      % saves a list of double occupations of all basis states
        %Index         % Basis index - useful for adressing Hamiltonian in many body basis
        Binary_Label    % label of binary representation
        Spin_Range_Data      % table of valid spin_range for particle number
        Categ

    end
    
    
    properties (Access = private)
        
    end
    
    methods
        function obj = Basis(N_site, Np_input, Sp_input, Docc_input, Num_of_Spin, sm_bi_sec)
            %N_site      ... number of system sites
            %Np_input    ... optional to restrict to certain particle numbers
            %Sp_input    ... optional to restrict to certain spin sectors (single spin +-1, thus total spin is integer numbers)
            %Docc_input  ... optional code to allow restrictions concerning double occupation / empty sites
            %Num_of_Spin ... optional to set up if spin is present (2) or not (1...spinless)
            %sm_bi_sec   ... convention in which order basis is set up from binary digits, standard is (u1,u2,u3|d1,d2,d3)
            %TODO ... this convention could be more general (u3,u2,u1,d3 ...; u4,d4,u3,d2, ...; d1,d2,d3,u1,..., 3 flavours
            %DONE ... estimate number of basis states
            if nargin < 1 || isempty(N_site); error('Input argument "Number of sites" is missing!'); end
            N_site = uint64(N_site);
           
            if nargin < 5 || isempty(Num_of_Spin), N_spin = uint64(2); else, N_spin = uint64(Num_of_Spin); end
            if nargin < 2 || isempty(Np_input), Np = 0:N_spin*N_site;
            else, Np_L = ismember(Np_input, 0:N_spin*N_site); 
                Np_input = uint64(Np_input);
                Np = Np_input(Np_L); 
            end %checks if particle input is in valid range (ignoring Docc)
            if nargin < 3 || isempty(Sp_input), Sp_all = true; else
                Sp_input = int64(Sp_input);
                Sp_all = false; 
            end
            if nargin < 4 || isempty(Docc_input); Docc = true; else, Docc = Docc_input; end
            if nargin < 6 || sm_bi_sec ~= 'u'; obj.Smaller_Binary_Sector_Used_For = 'down spin'; else
                obj.Smaller_Binary_Sector_Used_For = 'up spin'; end
            
            % Table for code of Docc_input
            %         |    0    |    1    |    2    |
            % ---------------------------------------
            %     1   |    y    |    y    |    y    |
            %     0   |    n    |    y    |    n    |
            %     2   |    y    |    y    |    n    |
            %     3   |    n    |    y    |    y    |
            %         |    y    |    n    |    y    | (equivalent to spinless Num_of_Spin = 1)
            %
            
            % optional - change order / conversion
            %   A - spin, site order of binary representation (horizontal order)
            %   b - order by spin / particle sectors (vertical order)
            %   in each segment the basis is sorted (useful for bisectional search)
            
            
            
            
            
            
            %NOTE: use int64 to describe basis states for up to 30 sites (60 electrons).
            %NOTE: use signed version since spin could be negative and negative numbers occur in the algorithm
            
            
            % logic predefined by given input
            % =============================================================
            % checking input (strict for Np_input, relaxed for Spin input (no warning or error))
            
            if Num_of_Spin == 1
                obj.Numel_Bas = sum(factorial(N_site)./(factorial(Np).*factorial(N_site-Np)));
            else
                
                if Docc == 0
                    if nargin >= 2 && ~isempty(Np_input)
                        warning('particle number input will be ignored since it has to be equal to number of sites since we want a Heisenberg basis');
                    end
                    Np = N_site;
                    
                    if ~Sp_all && all(~ismember(uint64(Sp_input),-N_site:2:N_site))
                        error('given spin sectors does not allow half filling')
                    end
                    
                    %calculating the number of basis elements
                    Up_range = 0:N_site;
                        if Sp_all == false
                            up_index = ismember((Sp_input+int64(Np))/2, int64(Up_range));
                            Up_range = uint64((Sp_input(up_index) + int64(Np))/2);
                            
                            
                        end
                    obj.Numel_Bas = sum(factorial(N_site)./(factorial(Up_range)  .* factorial( N_site - Up_range)));
                    
                elseif Docc == 1
                    % correcting Np not necessary -
                    % TODO: maybe implement Np correction totally here
                    
                    %calculating the number of basis elements - no spin check
                    numel_bas = 0;
                    for N_el = Np
                        Up_range = max(0,N_el-N_site):min(N_el, N_site);
                        if Sp_all == false
                            up_index = ismember((Sp_input+int64(N_el))/2, int64(Up_range));
                            Up_range = uint64((Sp_input(up_index) + int64(N_el))/2);
                            
                            
                        end
                        %comment: use nchoosek since rounding errors occur for very high numbers!!!
                        for ind_up = 1:numel(Up_range)
                            numel_bas = numel_bas + nchoosek(N_site, Up_range(ind_up)) * nchoosek(N_site, N_el-Up_range(ind_up));
                        end
                        %numel_bas = numel_bas + sum(factorial(N_site)^2./(factorial(Up_range) .* factorial(N_site - Up_range) .* factorial(N_el - Up_range) .* factorial(N_site - N_el + Up_range)));
                        
                    end
                    
                    obj.Numel_Bas = numel_bas;
                    
                elseif Docc == 2
                    if nargin >= 2 && ~isempty(Np_input)
                        if all(Np>N_site)
                            error('particle number has to be lower equal than number of sites (no double occupation)')
                        elseif any(Np> N_site)
                            warning('particle number has to be lower equal than number of sites (no double occupation)')
                            Np = Np(Np<=N_site);
                        end
                    else
                        Np = uint64(0:N_site);
                    end
                    
                    %calculating the number of basis elements - no spin check
                    numel_bas = 0;
                    for N_el = Np
                        Up_range = 0:N_el;
                        if Sp_all == false
                            up_index = ismember((Sp_input+int64(N_el))/2, int64(Up_range));
                            Up_range = uint64((Sp_input(up_index) + int64(N_el))/2);
                        end
                        numel_bas = numel_bas + sum(factorial(N_site)./(factorial(Up_range).*factorial(N_el - Up_range) .* factorial(N_site - N_el)));
                    end
                    obj.Numel_Bas = numel_bas;
                    
                    
                elseif Docc == 3
                    if all(Np<N_site)
                        error('particle number has to be greater equal than the number of sites (no vacant site)')
                    elseif any(Np<N_site)
                        warning('particle number has to be greater equal than the number of sites (no vacant site)')
                        Np = Np(Np>=N_site);
                    end
                    %calculating the number of basis elements
                    numel_bas = 0;
                    for N_el = Np
                        Up_range = (N_el - N_site) : N_site;
                        if Sp_all == false
                            up_index = ismember((Sp_input+int64(N_el))/2, int64(Up_range));
                            Up_range = uint64((Sp_input(up_index) + int64(N_el))/2);
                        end
                        numel_bas = numel_bas + sum(factorial(N_site)./(factorial(N_site - Up_range).*factorial(N_site - N_el + Up_range) .* factorial(N_el - N_site)));
                    end
                    obj.Numel_Bas = numel_bas;
                    
                end
            end
            
            
            % preallocation
            disp(['Number of basis elements: ',  num2str(obj.Numel_Bas)])
            
            bas_dez =   uint64(zeros(obj.Numel_Bas,1));
            n_part =     int8(zeros(obj.Numel_Bas,1));
            spin =       int8(zeros(obj.Numel_Bas,1));
            double_occ = int8(zeros(obj.Numel_Bas,1));
            
            count = 0;
            
            for N_el = Np
                Up_range = max(0,N_el-N_site):min(N_el, N_site);
                if Sp_all == false
                    up_index = ismember((Sp_input+int64(N_el))/2, int64(Up_range));
                    Up_range = uint64((Sp_input(up_index) + int64(N_el))/2);
                    
                    
                end
                
                
                if N_spin == 2
                    %Spin_range = 2*Up_range - N_el;
                    for N_up = Up_range
                        N_down = N_el - N_up;
                        
                        if N_site == 1
                            UP_site = uint64(zeros(1,N_up));
                        else
                            UP_site = nchoosek(0:N_site-1, N_up);
                        end
                        UP = sort(uint64(sum(2.^UP_site,2)));
                        
                        if Docc == 0
                            DOWN = 2^N_site -1 - UP;
                        else
                            if N_site == 1
                                DOWN_site = uint64(zeros(1,N_down));
                            else
                                DOWN_site = nchoosek(0:N_site-1, N_down);
                            end
                            DOWN = sort(uint64(sum(2.^DOWN_site,2)));
                        end
                        %convention to put UP spin to the left (in terms of binary exponentials)
                        %smaller binary sector is on the right side
                        if obj.Smaller_Binary_Sector_Used_For(1) == 'u'
                            DOWN = DOWN*2^N_site;
                        else
                            UP = UP * 2^N_site;
                        end
                        
                        if Docc == 0
                            basis_temp = UP(:) + DOWN(:);
                            
                            Double_Occup = uint64(zeros(size(basis_temp)));
                        else
                            [up,down] = meshgrid(UP,DOWN);
                            basis_temp = up(:)+down(:);
                            
                            bin = obj.dez2bin(basis_temp,2*N_site);
                            Double_Occup = sum(bin(:,1:N_site) & bin(:,N_site+1:end),2);
                            %Double_Occup = zeros(size(basis_temp));
                        end
                        
                        
                        %zero_occup missing
                        
                        %check special cases (see table above)
                        if  Docc == 2
                            L_del = Double_Occup > 0; %any site double occupied
                        elseif Docc == 3
                            L_del = sum(~bin(:,1:N_site) & ~bin(:,N_site+1:end),2) > 0; %any site not occupied
                        elseif Docc == 1 || Docc == 0
                            L_del = false(size(basis_temp));
                        end
                        %delete sites
                        basis_temp = basis_temp(~L_del);
                        Double_Occup = Double_Occup(~L_del);
                        
                        nnumel = numel(basis_temp);
                        
                        bas_dez(count + (1:nnumel)) = basis_temp;
                        double_occ(count + (1:nnumel)) = Double_Occup;
                        spin(count + (1:nnumel)) = int8(ones(size(basis_temp))) * (2*int8(N_up) - int8(N_el));
                        n_part(count + (1:nnumel)) = uint64(ones(size(basis_temp))) * N_el;
                        
                        count = count + nnumel;
                        %                         elseif Docc == 0
                        %                             UP_site = nchoosek(0:N_site-1, N_up);
                        %                             UP = sum(2.^UP_site,2);  %same
                        %                               %new
                        %
                        %                             %convention to put UP spin to the left (in terms of binary exponentials)
                        %                             UP = UP * 2^N_site;  %same
                        %                             obj.Bas_Dez = [obj.Bas_Dez; UP(:)+DOWN(:)]; %quite the same
                        %                             obj.Double_Occ = [obj.Double_Occ;zeros(size(UP(:)))];  %easier
                        %                             obj.spin = [obj.spin; ones(size(up(:)))*(2*N_up - N_el)]; %same
                        %                             obj.N_part = [obj.N_part; ones(size(up(:)))*N_el];  %same
                        %                         end
                        
                    end
                    
                    if strcmp(obj.Smaller_Binary_Sector_Used_For, 'up spin')
                        obj.Binary_Label = [strcat( repmat({'d'},1,N_site),cellstr(num2str((1:N_site)'))'), ...
                            strcat(repmat({'u'},1,N_site), cellstr(num2str((1:N_site)'))')];
                    else
                        obj.Binary_Label = [ strcat( repmat({'u'},1,N_site), cellstr(num2str((1:N_site)'))'), ...
                            strcat( repmat({'d'},1,N_site), cellstr(num2str((1:N_site)'))')];
                    end
                elseif N_spin == 1
                    % feature of spinless fermions
                    if N_site == 1
                        part_exp = uint64(zeros(1,N_el));
                    else
                        part_exp = nchoosek(0:N_site-1, N_el);
                    end
                    part = sum(2.^part_exp,2);
                    
                    nnumel = numel(part(:));
                    
                    bas_dez(count + (1:nnumel)) = part(:);
                    n_part(count + (1:nnumel)) = uint64(ones(size(part)))*N_el;
                    %double_occ(count + (1:nnumel)) = zeros(size(part)); % not existent - already with zeros initialized
                    %spin(count + (1:nnumel)) =  zeros(size(part)); % not existent - already with zeros initialized
                  
                    count = count + nnumel;
                    
                    obj.Binary_Label = strcat( repmat({'s'},1,N_site), cellstr(num2str((1:N_site)'))'); %num2cell(1:N_site);
                end
                
            end
            if numel(bas_dez) == 0
                error('no basis set found for the given parameters')
            end
            
            obj.Categ = Category({'Index', 'I'}, uint64(1:obj.Numel_Bas).');
            obj.N_site = double(N_site);
            obj.N_spin = double(N_spin);
            
            
            %initalize Category
            obj.Categ = obj.Categ.f_append({'Bas_Dez', 'B'}, bas_dez);
            obj.Categ = obj.Categ.f_append({'N_part', 'N'}, n_part);
            
            if obj.N_spin > 1
                obj.Categ = obj.Categ.f_append({'Spin', 'S'}, spin);
                obj = obj.f_generate_spin_range;
            end
            obj.Categ = obj.Categ.f_append({'Double_Occupancy', 'D'}, double_occ);
            
            
            obj.Numel_Bas = double(obj.Numel_Bas);
            
            
         
        end
        
        % ================================================================
        % ------------------ getter functions!!!! ------------------------
        function index = Index(obj, varargin)
            index =  obj.Categ.f_access('Index', varargin{:});
        end
        
        function bas_dez = Bas_Dez(obj, varargin)
            bas_dez =  obj.Categ.f_access('Bas_Dez', varargin{:});
        end
        
        function n_part = N_part(obj, varargin)
            n_part =  obj.Categ.f_access('N_part', varargin{:});
        end
        
        function spin = Spin(obj, varargin)
            if obj.N_spin == 1
                error('spinless system!')
            end
            spin = obj.Categ.f_access('Spin', varargin{:});
        end
        
        function double_occ = Double_Occ(obj, varargin)
            double_occ = obj.Categ.f_access('Double_Occupancy', varargin{:});
        end
        
        % ===============================================================
        
        % ===============================================================
        % -------------------address site representation ----------------------------
        
        
        function bas_bin = f_bin(obj, varargin)
            %NEW function returns binary representation of spin, uses categories to filter
            %NOTE: spin_type feature (lower or upper spins) not used!
            if nargin == 1
                bas_dez = obj.Bas_Dez;
            else
                bas_dez = obj.Bas_Dez(varargin{:});
            end
            
            bas_bin = obj.dez2bin(bas_dez, obj.N_spin * obj.N_site);
            
        end
        
        % ===============================================================
        % -------------------print functions ----------------------------
        
        function [] = s_print_full(obj)
            %returns table of basis, filter options: not implemented yet
            if obj.N_spin == 1
                disp(table(obj.Bas_Dez, obj.N_part, obj.Double_Occ, 'variablenames',{'Bas_dez', 'N', 'Double_Occ'}))
            else
                disp(table(obj.Bas_Dez, obj.N_part, obj.Spin, obj.Double_Occ, 'variablenames',{'Bas_dez', 'N', 'Spin', 'Double_Occ'}))
            end
        end
        
        function [Tab] = s_print(obj)
            if obj.N_spin == 1
                [A, ~, freq] = unique([obj.N_part], 'rows');
                numb = histc(freq,1:size(A,1));
                Tab =   array2table([A,numb], 'variablenames',{'N_part', 'N_Basis_Elements'});
            else
                [A, ~, freq] = unique([obj.N_part, obj.Spin], 'rows');
                numb = histc(freq,1:size(A,1));
                Tab =   array2table([A,numb], 'variablenames',{'N_part', 'Spin', 'N_Basis_Elements'});
            end
            disp(Tab)
            
        end
        
        
        % =================================================================
        % ---------------------spin range functionality ------------------
        function spin_rg = func_spin_range(obj, N_el)
            %calculates the spin range for a single N_el input
            
            if obj.N_spin == 2
                spin_rg = 2*(max(0,N_el-obj.N_site):min(N_el,obj.N_site)) - N_el;
            else
                warning('another N_spin as two defined, N_spin == 1 shouldnot reach this point.')
                spin_rg = 0;
            end
            
        end
        
        function obj = f_generate_spin_range(obj)
            %generates all spin_ranges (only called when N_spin ~= 1)
            for N_el = 0:(obj.N_site*obj.N_spin)
                obj.Spin_Range_Data{N_el+1} = obj.func_spin_range(N_el);
            end
        end
        
        function spin_range = f_spin_range(obj, N_el, Spin)
            %three features:
            % - returns spin range from stored cell (since recomputation is too expensive)
            %   cell is produced by f_generate_spin_range 
            %       - depending on particle sector (N_el)
            %       - for all particles
            % - checks if Spin is in spin_range
            
            % NOTE: works for scalar input only at the moment
            %
            if obj.N_spin == 1
                error('spinless system!')
                %spin_range = 0; % NOTE: old version
            end
            if nargin < 3 || isempty(Spin)
                
                if nargin < 2 || isempty(N_el)
                    % return all possible spins for all particle sectors
                    %spin_range = obj.Spin_Range_Data{}
                    spin_range = unique(cell2mat(obj.Spin_Range_Data));
                else
                    if N_el < 0 || N_el > obj.N_spin*obj.N_site
                        spin_range = [];
                    else
                        spin_range = obj.Spin_Range_Data{N_el+1};
                    end
                end
               
            else
                %NOTE: the check part works also for an input vector Spin
                spin_range = all(ismember(Spin, obj.Spin_Range_Data{N_el+1}));
            end

            
        end
        
        % ======================== symmetric operations ===================
        function index = f_symmetric_action(obj, operation, varargin)
            bas_bin = obj.f_bin(varargin{:});
            bas_ref = obj.Bas_Dez(varargin{:});
            
            if iscell(operation)
    
                bas_new = bas_bin(:,[operation{1}, obj.N_site +operation{2}]);
            else
                if obj.N_spin == 1
                    bas_new = bas_bin(:,operation);
                else
                    bas_new = bas_bin(:,[operation, obj.N_site + operation]);
                end
            end
            
            [~,index] = ismember(obj.bin2dez(bas_new), bas_ref);

            
        end
   
    
        % =================================================================
        % ---------------------not used functions? - check ----------------
        
        function bin_index = f_bin_index(obj, c_site, spin)
            % returns the index of the binary Basis corresponding to desired site and spin
            % returns 0 if input is not valid
            %TODO: check if this works also for zero spin
            % check where this function is needed!!!
            warning('not sure where this function is used!')
            Lerr = c_site <= 0 | c_site > obj.N_site;
            Lerr2 = ~ismember(spin,[-1,1]);
            
            L1 = obj.Smaller_Binary_Sector_Used_For(1) == 'd';
            L2 = spin == 1;
            
            
            shift = xor(L1,L2)*obj.N_site;
            bin_index = c_site + shift;
            bin_index(Lerr | Lerr2) = 0;
            
        end
        
        
        function bas_bin = g_bin_by_index(obj, index)
            warning('not sure where this function is used!')
            bas_bin = obj.dez2bin(obj.Bas_Dez(index), obj.N_site*obj.N_spin);
        end
        
        function bas_bin = g_bin_by_dez(obj, dez)
            warning('not sure where this function is used!')
            bas_bin = obj.dez2bin(dez, obj.N_site*obj.N_spin);
        end
        
        % =================================================================
        % ------------------old functions ---------------------------------
%         function [bas, N, spin, Double_Occ, index] = f_N(obj,Np_input,Sp_input,Docc_input)
%             error('not used any more!!! - use Basis.Bas_Dez instead')
%             % Particle sectorBasis(N_site, Np_input, Sp_input, Docc_input, Num_of_Spin, sm_bi_sec), Spin sector, Double occupation
%             
%             if nargin < 2 || isempty(Np_input), L_N = true;
%             else L_N = ismember(obj.N_part, Np_input); end  %if ismember is too slow use internal function ismembc (in combination with sort)
%             if nargin < 3 || isempty(Sp_input), L_S = true;
%             else L_S = ismember(obj.Spin,Sp_input); end
%             if nargin < 4 || isempty(Docc_input), L_D = true;
%             else L_D = ismember(obj.Double_Occ,Docc_input); end
%             L = L_N & L_S & L_D & true(size(obj.Bas_Dez));
%             
%             bas = obj.Bas_Dez(L);
%             N = obj.N_part(L);
%             spin = obj.Spin(L);
%             Double_Occ = obj.Double_Occ(L);
%             index = obj.Index(L);
%         end
%         
%         function index = f_ind(obj, Np_input, Sp_input, Docc_input)
%             % returns basis index, filter options: part_num, spin, Double_Occ
%             error('not used any more - use Basis.Index instead!')
%             if nargin < 2 || isempty(Np_input), L_N = true;
%             else L_N = ismember(obj.N_part, Np_input); end  %is ismember is too slow use internal function ismembc (in combination with sort)
%             if nargin < 3 || isempty(Sp_input), L_S = true;
%             else L_S = ismember(obj.Spin,Sp_input); end
%             if nargin < 4 || isempty(Docc_input), L_D = true;
%             else L_D = ismember(obj.Double_Occ,Docc_input); end
%             L = L_N & L_S & L_D & true(size(obj.Bas_Dez));
%             
%             index = obj.Index(L);
%         end
%         function bas_bin = f_bin_(obj, Np_input, Sp_input, Docc_input, spin_type)
%             %returns binary representation of basis, filter options: part_num, spin, doubl_occ, spin_type
%             error('not used any more - use Basis.f_bin instead!')
%             if nargin < 2 || isempty(Np_input), Np_input = 0:obj.N_spin * obj.N_site; end
%             if nargin < 3 || isempty(Sp_input), Sp_input = []; end
%             if nargin < 4 || isempty(Docc_input), Docc_input = []; end
%             if nargin < 5 || isempty(spin_type) || (lower(spin_type(1)) ~= 'u' && lower(spin_type(1)) ~= 'd') || obj.N_spin ~=2,
%                 spin_type = 'b'; end
%             
%             bas_bin = obj.dez2bin(obj.f_N(Np_input, Sp_input, Docc_input), obj.N_spin * obj.N_site);
%             if spin_type == 'd'
%                 bas_bin = bas_bin(:,1:obj.N_site);
%             elseif spin_type == 'u'
%                 bas_bin = bas_bin(:,1:obj.N_site);
%             end
%             
%         end
        
    end
    
    
    
    methods(Static, Access = private)
        function base_bin = dez2bin(base,Length)
            base_bin = uint64(zeros(numel(base),Length));
            for i = 1:Length
                I = base>=2^(Length-i);
                base_bin(I,i) = 1;
                base(I) = base(I)-2^(Length-i);
                
            end
        end
        
        function base = bin2dez(base_bin)
            base = bin2dec(num2str(base_bin));
        end
    end
    
    
    
    
    
    
end