function [K,I,Sig, Bat, Tab_cut, status, column_sum, ev_min, choi, Sig_result] = born_markov_full_dynamic_old(En, contact_site, Tab, Mu, Coupling, Beta, Bat, condor, Energy_cut, Energy_cut_2, run_id, min_en_ind, eigenvalue_tolerance)

if ismember('S', Tab.Order); spin_order = 'S'; else, spin_order = []; end
%Qtab = Table(En, ['N', spin_order]);
if nargin < 12 || isempty(min_en_ind), min_en_ind = min(50, sum(energ_sort < (Energy_cut + min(energ_sort))/2)); end

if nargin < 11 || isempty(run_id), run_id = ''; end
if nargin < 10 || isempty(Energy_cut_2), Energy_cut_2 = inf; end
if nargin < 9 || isempty(Energy_cut), Energy_cut = []; end
% L = Qtab.N_part <= 3;
% Qtab.N_part = Qtab.N_part(L);
% Qtab.Spin = zeros(sum(L),0);
% Qtab.Energies = zeros(sum(L),0);
% Qtab.Block_Size = Qtab.Block_Size(L);
% Qtab.Index = Qtab.Index(L);
% Qtab.Vector_Index = Qtab.Vector_Index(L);
% Qtab.Diag_Index = Qtab.Diag_Index(L);
% for k = 1:3
% Qtab.Category(k).Value = Qtab.Category(k).Value(L);
% end

%TODO: define as input variable: 
sigma_tolerance = 10^-10;

if ismac
    % Code to run on Mac plaform
    SLASH = '/';
elseif isunix
    % Code to run on Linux plaform
    SLASH = '/';
elseif ispc
    % Code to run on Windows platform
    SLASH = '\';
else
    disp('Platform not supported')
end



mes_flag = condor.mes_flag;
if ~isfield(condor, 'setup_id')
    condor.setup_id = '';
end
if ~isfield(condor, 'save_intermediate_K')
    condor.save_intermediate_K = false;
end

Tab_cut = Tab;

%TODO: change input for all files (maybe make a conversion file, so that old code still works
% gam_vec(baths, spin) -> Coupling{baths, spin}
% mu_vec(left,right) -> Mu(baths, spin)
% beta(1) -> Beta(baths, spin)
% Bat -> new Bath class (think of it)

%conversion
if numel(Beta) == 1 && (iscell(Coupling) == false)
    warning('conversion of gam_vec, beta and mu')
    N_site = En.Basis.N_site;
    [N_bath_inp, N_spin_inp] = size(Coupling);
    Coupling_big = zeros(N_site*N_spin_inp,N_bath_inp);
    for ii = 1:N_bath_inp
        for jj = 1:N_spin_inp
            Coupling_big((ii-1)*N_site + contact_site(ii),  jj) = Coupling(ii,jj);
        end
    end
    Coupling = mat2cell(Coupling_big,N_site*ones(1,N_spin_inp) ,ones(N_bath_inp,1));
    %Syntax for coupling: Coupling{Bath, Spin}(System Site, Bath orbital) !!!!
    Beta = ones(size(Coupling))*Beta;
    Mu = repmat(Mu(:), 1, size(Coupling,2));
end

N_bath = size(Coupling,1);
if (Bat.Numel_Baths ~= N_bath) || (size(Mu,1) ~= N_bath) || size(Beta,1) ~= N_bath, error('number of baths is not consistent'); end

rel_site = 1:En.Basis.N_site; % to which site there is coupling: relevant sites
rel_site = rel_site(ismember(rel_site, contact_site)); %each counted only once

% loop over all blocks

% lookup table to have an indexing scheme which elements have to be taken
if ismember('E', Tab.Order) %search in sectors of Q-matrices
    % the other cases are, that there is no spin, spin not conserved, then 
    % Q matrices don't have a spin!
    %index_translate = false(1,Qtab.Numel_Vector_Elements);
    %index_translate(Tab.f_index_translation(1:Tab.Numel_Vector_Elements, Qtab)) = true;
    %sel_block = Qtab.f_block(index_translate);
    secular_approximation = true;
else % ismember('A', Tab.Order)
    %TODO: find better way to count number of blocks
    %index_translate = true(1,Tab.Numel_Vector_Elements);
    %sel_block = Tab.f_block(index_translate);
    secular_approximation = false;
end

if En.Basis.N_spin == 1
    spinless = true;
    c_spin_range = -1;
else 
    spinless = false;
    c_spin_range = [-1,1];
end

if ismember('S', Tab.Order)
    spin_conserved = true;
    spin_tag = 'S';
else
    spin_conserved = false;
    spin = [];
   
    spin_tag = 'S';
end

%flag to identify status of calculation
status.no_solution = 0; %indicates if no solution (=1) is found
status.multiple_solutions = 0; % indicates number of solutions (zero eigenvalues)
status.not_same_argument = 0; %indicates whether blocks have different complex argument
status.current_missmatch = 0; % indicates whether current in and out are conserved for the ?? method
status.column_sum = 0; % indicates the column sum
status.positive_eigenvalues = 0; % indicates
status.eigenvalue_tolerance = eigenvalue_tolerance;
status.minimal_value = inf;
status.energy_cut = inf;
status.error_text = '';
status.choi_herm = inf;
status.choi_neg_ew_ratio = 0;
status.choi_mean_ew = -inf;

symmetrized_version = false;
include_sigma_energy = false;

 [column_sum, ev_min, choi] = deal(0);

%TODO: maybe reduce sparse size due to energy cut
%K = sparse(sum(Tab.Block_Size.^2), sum(Tab.Block_Size.^2));

%blocks that are used (energy cut)
if ismember('E', Tab.Order)
    energy_tag = 'E';
    table_tag = ['NS', 'E'];
elseif ismember('A', Tab.Order)
    energy_tag = 'A';
    table_tag = ['NS', 'A'];
end


if ~isempty(Energy_cut)
    energ_index = Tab.Block_Index(table_tag, [], [], {-1000000, Energy_cut});
else
    energ_index = Tab.Block_Index;
end

[I1, I2, W] = deal(zeros(sum(Tab.Block_Size(energ_index).^2)*En.Basis.N_spin*En.Basis.N_site*4,1));
sparse_index = 0;

if mes_flag
    disp([num2str(numel(energ_index)), ' blocks of ', num2str(numel(Tab.Block_Index)) ' used'])
end
sigma_index = [];

count = 0;

% sort the blocks by their energy
if ismember('A', Tab.Order)
    [energ_sort, energ_sort_index] = sort(Tab.Average_Energies);
else
    [energ_sort, energ_sort_index] = sort(Tab.Energies);
end
%[energ_sort, energ_sort_index]
energ_sort_index(energ_sort  > Energy_cut) = []; %delete blocks that exceed Energy_cut

final_block_index = numel(energ_sort_index);
break_set = final_block_index; %[min_en_ind:10:numel(energ_sort_index), final_block_index];
col_sum = zeros(size(break_set));
col_sum_index = 0;

X = cell( numel(energ_sort_index) , N_bath, numel(c_spin_range));
X_sym = cell( numel(energ_sort_index) , N_bath, numel(c_spin_range));
Bath3D = cell(N_bath,1);
Bath3D_2 = cell(N_bath,1);

run_K_flag = true;

for n_block_index = 1:numel(energ_sort_index) 
    if run_K_flag
        %for n_block_index = 5
        n_block = energ_sort_index(n_block_index);
        if mes_flag
            if any(ismember('E', Tab.Order))
                energ_tag = [', Energy: ', num2str(Tab.Energies(n_block))];
            elseif any(ismember('A', Tab.Order))
                energ_tag = [', Energy: ', num2str(Tab.Average_Energies(n_block))];
            else
                energ_tag = [];
            end
            
            
            disp(['Block number: ' num2str(n_block), '/', num2str(size(Tab.N_part,1)), ', Block: ', num2str(n_block_index), '/', num2str(numel(energ_index)) , ', Blocksize: ', num2str(Tab.Block_Size(n_block)), ...
                energ_tag])
            %disp(['Sparse index: ', num2str(max(sparse_index))])
        end
        
        if ismember('N', En.Structure)
            n_part = Tab.N_part(n_block);
        else
            n_part = [];
        end
        if spin_conserved %spin is preset to [] (not spin preserved and spinless case
            spin = Tab.Spin(n_block);

        end
        
        index1 = Tab.Vector_Index(n_block);
        
        block_size = Tab.Block_Size(n_block);
        X1 = zeros(block_size.^2);
        
        index1_trans = reshape(1:block_size^2, block_size,block_size)';
        %X1 = zeros(numel(index1));
        % TODO: optimization for secular approximation!!!
        %for secular approximation and dealing with a density operator
        % X1 can be smaller:
        % X1 = zeros(Tab.Block_Size(n_block));
        % perform later the Kronecker tensor product:
        % K(index1, index1) = K(index1, index1) + kron(eye(Tab.Block_Size(n_block)),X1)
        
        
        %block_index needed just for secular_approximation
        %TODO: replace by new Tab index functions!!!!
        %n_block_first = find(Tab.N_part == n_part & L_S, 1, 'first'); %TODO: check if this works
        %block_index = sum(Tab.Block_Size(n_block_first:(n_block-1))) + (1:Tab.Block_Size(n_block));
        %Er = En.g_En(n_part, spin);
        %Er = Er(block_index); %should be one energy if secular approximation is on
        
        % NOTE: more efficient (check if also fast)
        list_index = Tab.List_Index('I', n_block);
        
        %think about if Average Energy would be better?
        Er = En.Energies('I', list_index);
        sub_block_index = En.Subindex('NS',  {n_part, spin}, 'I', {list_index});
        
        %index of entries for elements in particle/spin block with certain energy (no need of Qtab here)
        %needed for secular approximation (work in energy blocks)
        
        %
        %NOTE: maybe use Average energy???
        %NOTE: relevant for Qmatrices
        
        %initialize X for current calculation
        %     for alpha = 1:N_bath
        %         for c_spin = c_spin_range
        %             X{n_block, alpha, c_spin/2+1.5} = 0;
        %         end
        %     end
        [X{n_block, :, :}] = deal(0);
        [X_sym{n_block, :, :}] = deal(0);
        for dag = -1:2:1
            
            for c_spin = c_spin_range %for spinless c_spin_range is just -1
                %%
                if spinless
                    c_spin = [];
                end
                
                count = count +1;
                
                %disp(['round: ', num2str(count)])
                %NOTE: querry for virtual excitations: a^s a^_s \rho  (Energy_cut_2) case A)
                %NOTE: for second part only Energy_cut relevant a^s \rho a^_s        case B)
                if numel(Tab.Block_Index(table_tag, n_part+dag, spin + dag*c_spin, {-inf, Energy_cut_2})) >= 1
                    
                    %TODO: Test with 'A' instead of E
                    energy_subindex = En.Subindex('NS', {n_part + dag, spin + dag*c_spin}, energy_tag, {{-inf, Energy_cut}});
                    
                    %energy_subindex = En.f_energy_subindex(n_part + dag, spin + dag*c_spin, {-inf, Energy_cut});
                    next_block_sizes = Tab.Block_Size(table_tag, n_part + dag, spin + dag*c_spin, {-inf, Energy_cut});
                    
                    %TODO: Test with 'A' instead of E
                    energy_subindex_virtual = En.Subindex('NS',{n_part + dag, spin + dag*c_spin}, energy_tag, {{-inf, Energy_cut_2}});
                    %energy_subindex_virtual = En.f_energy_subindex(n_part + dag, spin + dag*c_spin, {-inf, Energy_cut_2});
                    
                    %TODO: Test with 'A' instead of E
                    El = En.Energies(['NS', energy_tag], n_part + dag, spin + dag*c_spin, {-inf, Energy_cut_2});
                    [Enr, Enl] = meshgrid( Er, El );
                    
                    A3 = zeros(block_size^2 , sum(next_block_sizes.^2) ); %relevant for case B)
                    
                    Enr2 = Enl.'; Enl2 = Enr.'; 
                    
                    %NOTE: could be constant for constant Gammas
                    
                    for alpha = 1:N_bath
                        mu_alpha = Mu(alpha, Bat.Spin_fun(c_spin));
                        beta_alpha = Beta(alpha, Bat.Spin_fun(c_spin));
%                         
%                         Enl = linspace(-2.1,2,200);
%                         Enr = zeros(size(Enl));
%                       
                        % first term, form: A*(C.*B)*sigma; A = c_\kappa^{\overline{s}}, B = c_\mu^s
                        if symmetrized_version
                            %TODO: check correct dag or -dag operator, maybe it is possible to exchange both
                            %without a significant change.
                              Bath3D{alpha} = Bat.f_gamma(-(Enl(:)-Enr(:)), mu_alpha, beta_alpha, alpha, c_spin, -dag);
                              Bath3D_sym{alpha} = Bat.f_gamma(Enr(:) - Enl(:), mu_alpha, beta_alpha, alpha, c_spin, dag);
                              precalc_warning = 0;
%                             figure
%                             plot(Enl, squeeze(real(C)))
%                             hold on
%                             plot(Enl, squeeze(imag(C)))
                        else
                            if false
                                %TODO: check formulas for different dagger, negative dagger sign !
                                C_real = Bat.f_gamma(-(Enl(:)-Enr(:)), mu_alpha, beta_alpha, alpha, c_spin, -dag);
                                
                                C_imag = Bat.f_sigma(-(Enl(:)-Enr(:)), mu_alpha, beta_alpha, alpha, c_spin, -dag);
                                C = 1/2*(C_real + C_imag);
                                precalc_warning = 0;
                                Bath3D{alpha} = C;
                            end
                            
                            [Bath3D{alpha}, Bat, precalc_warning] = Bat.f_lookup(-(Enl(:)-Enr(:)), mu_alpha, beta_alpha, alpha, c_spin, dag);
%                             E = squeeze(Bat.f_calculate_detail(-(Enl(:)-Enr(:)), mu_alpha, beta_alpha, alpha, c_spin, dag));
%                             %F = squeeze(Bat.f_flat_band_calculate(-(Enl(:)-Enr(:)), mu_alpha, beta_alpha, alpha, c_spin, dag));
% 
% 
%                             D = squeeze(Bath3D{alpha})
%                             plot(Enl, real(E), Enl, imag(E))
                        end
                        
                        if precalc_warning
                            
                            disp(['N: ', num2str(n_part), ', S: ', num2str(spin), ', N": ', num2str(n_part + dag), ', S": ', num2str(spin + dag*c_spin)])
                            %disp(num2str(-Enl(1) - Enr(:)))
                        end
                               
                        
                        %TODO: check if coupling_mu coupling_kappa sequence is correct
                        %Bath3D = Bat.f_calculate(-(Enl2(:)-Enr2(:)), mu_alpha, beta_alpha, alpha, c_spin, -dag);
                        
                        if symmetrized_version
                            Bath3D_2{alpha} = Bat.f_gamma(-(Enl2(:) - Enr2(:)), mu_alpha, beta_alpha, alpha, c_spin, dag);
                            precalc_warning = 0;
                            
                            
                        else
                            if false
                                C_real = Bat.f_gamma(-(Enl2(:) - Enr2(:)), mu_alpha, beta_alpha, alpha, c_spin, +dag);
                                
                                C_imag = Bat.f_sigma(-(Enl2(:) - Enr2(:)), mu_alpha, beta_alpha, alpha, c_spin, +dag);
                                C = 1/2*(C_real + C_imag);
                                precalc_warning = 0;
                                Bath3D_2{alpha} = C;
                            end
                            
                            
                            [Bath3D_2{alpha}, Bat, precalc_warning] = Bat.f_lookup(-(Enl2(:)-Enr2(:)), mu_alpha, beta_alpha, alpha, c_spin, -dag);
                        end
                        
                        if precalc_warning
                            disp('precalc warning')
                        end
                        
                    end
                    
                    
                    
                    
                    for site_index_mu = 1:numel(rel_site)
                        c_site_mu = rel_site(site_index_mu);
                        
                       
                        
                        B = En.g_Qmat(dag, n_part, spin, c_site_mu, c_spin, Energy_cut_2);
                        B_virt = B(energy_subindex_virtual,sub_block_index);
                        
                        for site_index_kappa = 1:numel(rel_site)
                            c_site_kappa = rel_site(site_index_kappa);
                            % TODO: check if there is coupling from mu to kappa, or if loop run is not necessary
                            
                            
                            
                            
                            
                            
                            %structure: C_
                            %associated parameters
                            % Q-matrix and eigenenergies of sectors
                            
                            %dag: C and the operator C is acting on (B), have the same s (dag)
                            % A_kappa^{-dag} *(C.*(B_mu^dag * sigma)
                            
                            %                         B = En.g_Qmat(dag, n_part, spin, c_site_mu, c_spin, Energy_cut);
                            %                         B = B(:,block_index);
                            %
                            %                         A = En.g_Qmat(-dag,n_part + dag,spin + dag*c_spin, c_site_kappa, c_spin, Energy_cut);
                            %                         A = A(block_index, :); % here opposite, since A is kind of adjoint
                            %
                            %
                            %                         % Energy block from the left side (in the middle)
                            %                         El = En.g_En(n_part + dag, spin + c_spin * dag);
                            %
                            
                            %introduce a stronger Energy cut:
                            
                            
                            
                            A = En.g_Qmat(-dag,n_part + dag,spin + dag*c_spin, c_site_kappa, c_spin, Energy_cut_2);
                            A_virt = A(sub_block_index, energy_subindex_virtual); % here opposite, since A is kind of adjoint
                            
                            
                            % Energy block from the left side (in the middle)
                            
                            
                            
                            
                            
                            C = zeros(size(Enl,1),size(Enl,2), N_bath);
                            C_sym = zeros(size(Enl,2), size(Enl,1), N_bath);
                            for alpha = 1:N_bath
                                coupling_kappa = Coupling{alpha,Bat.Spin_fun(c_spin)}(c_site_kappa,:);
                                coupling_mu = Coupling{alpha,Bat.Spin_fun(c_spin)}(c_site_mu,:);
                                
                                
                                
                                % alpha for Bat means different bath correlation function for different baths
                                %TODO: optional: fasten up if different bath correlation functions are the same for
                                %different baths - same is valid for spin (would be better than copying whole Bat class
                                %four times - at the moment this is convenient! - check could also happen in Bat
                                %function of in rescaling index (if condition on Bath class Property (spin indep, bath
                                %symmetric bath)
                                
                                %TODO: check if coupling_mu coupling_kappa sequence is correct
                                
                                
                                %Bath3D = Bat.f_calculate(-(Enl(:)-Enr(:)), mu_alpha, beta_alpha, alpha, c_spin, dag);
                                %this is a 3 dimensional array (orbitals i and j for first dimensions), \Delta Energy
                                %for 3rd dimension. Aim is to perform multiple matrix multiplication for each energy
                                %and then to retrieve energy matrix (as done before)
                                
                                % TODO: Think about dimensions of energy!!!, Dimensions of Bath3D (
                                %C = C + reshape(squeeze(mmat(mmat(coupling_kappa,Bath3D{alpha}), coupling_mu')), size(Enl,1),size(Enl,2));
                                
                                C(:,:,alpha) = reshape(squeeze(mmat(mmat(coupling_kappa,Bath3D{alpha}), coupling_mu')), size(Enl,1),size(Enl,2));
                                X{n_block, alpha, Bat.Spin_fun(c_spin)} =X{n_block, alpha, Bat.Spin_fun(c_spin)} + dag * A_virt*(C(:,:,alpha).*B_virt);
                                %TODO: if sigma is not included X can be used to calculate K faster!!! -> implement
                                if symmetrized_version
                                    C_sym(:,:,alpha) = reshape(squeeze(mmat(mmat(coupling_kappa,Bath3D_sym{alpha}), coupling_mu')), size(Enl,2),size(Enl,1));
                                    X_sym{n_block, alpha, Bat.Spin(fun(c_spin))} = X_sym{n_block, alpha, Bat.Spin_fun(c_spin)} -1i/4*(A_virt*(C(:,:,alpha).*B_virt) + (C_sym(:,:,alpha).*A_virt)*B_virt);
                                end
                            end
                            
                            if any(isnan(C))
                                error('there should not be nans')
                            end
                            
                            %-------------------first term-------------------
                            % form: A*(C.*B)*sigma;
                            % C = C^s, B = B^s; A = A^-s
                            
                            % TODO: simplification for sigma quadratic and secular approximation possible:
                            % X1 = A*(C.*B);
                            
                            [~,sig2] = deal(Tab.Block_Size(n_block)); %should be squared
                            
                            [~,sizB2] = size(B_virt);
                            C = sum(C,3);
                            C_sym = sum(C_sym,3);
                            
                            BB = kron(eye(sig2), B_virt);
                            
                            if include_sigma_energy
                                CC1 = C(:);
                                CC = repmat(CC1, 1, sizB2 * sig2);
                            else
                                CC = kron(eye(sig2), C);
                                CC_sym = kron(eye(sig2),C_sym);
                            end
                            AA = kron(eye(sig2), A_virt);
                            
                            %disp([num2str(n_part),', S: ', num2str(spin),', dag: ',  num2str(dag), ', c_spin: ' , num2str(c_spin)])
                            %disp(['N: ', num2str(n_part), ' S: ', num2str(spin), ', dag: ', num2str(dag), ...
                            %    ', c_spin: ', num2str(c_spin), ', mu: ', num2str(site_index_mu), ', kappa: ', num2str(site_index_kappa)])
                            
                            %disp(AA*(CC.*BB))
                            if symmetrized_version
                                X1 = X1 + 1/2*( AA*(CC.*BB) + (AA.*CC_sym)*BB);
                            else
                                X1 = X1 + AA*(CC.*BB); % diagonal block element
                            end
                            %AA*(CC.*BB)
                            % TODO: think of complex conjugate
                            
                            
                            
                            % --------------------second term-----------------
                            
                            if numel(energy_subindex) >= 1
                                %TODO: reduce energies which are in between Energy_cut and Energy_cut_2
                                
                                %            off diagonal block element (index1, index2)
                                %  the question here is, what's more expensive: (a) second C_func or (b) second En.g_Qmat + En.En query
                                % (b): to recycle the Q-matrix we change dag -> -dag, energies are also the same but flipped:
                                % structure: C^-s .* A' * sigma * A
                                % sigma in N+dag sector
                                Enl2 = Enr.';
                                %B = En.g_Qmat(dag, n_part - dag, spin - dag*c_spin, rel_site, c_spin); % option (a) to
                                %recylce C_func (if there is an integration over the spectral function of the bath)
                                
                                
                                %C2 =  Bat.f_interp(-(Enl2-Enr2), mu_alpha, beta, -dag, Gamma);
                                C2 = zeros(size(Enl2));
                                for alpha = 1:N_bath
                                    
                                    coupling_kappa = Coupling{alpha,Bat.Spin_fun(c_spin)}(c_site_kappa,:);
                                    coupling_mu = Coupling{alpha,Bat.Spin_fun(c_spin)}(c_site_mu,:);
                                    
                                    
                                    %Bath3D = Bat{alpha,c_spin}.f_lookup(-(Enl2(:)-Enr2(:)), mu_alpha, beta_alpha, -dag);
                                    C2 = C2 + reshape(squeeze(mmat(mmat(coupling_kappa,Bath3D_2{alpha}), coupling_mu')), size(Enl2,1),size(Enl2,2));
                                    
                                end
                                if any(isnan(C2))
                                    error('there should not be nans')
                                end
                                
                                
                                
                                
                                %-------------------second term-------------------
                                % form: (C.*(B*sigma))*A; or (C.*B)*sigma * A
                                % C = C^s, A = A^s; A' = A^-s
                                
                                %get next block indices (due to energy cut you don't get all!)
                                %NOTE: taking full block and cutting then is equal to cut beforehand
                                
                                next_block_ind = Tab.Block_Index(table_tag, n_part + dag, spin + dag*c_spin, {-inf, Energy_cut});
                                B2 = A(sub_block_index, energy_subindex);
                                A2 = B(energy_subindex,sub_block_index);
                                
                                AI = cell(1,numel(next_block_ind));
                                blk= cell(1,numel(next_block_ind));
                                intra_block_size = Tab.Block_Size(next_block_ind);
                                for k = 1:numel(next_block_ind)
                                    blk_temp = (1:intra_block_size(k)) + sum(intra_block_size(1:k))-intra_block_size(k);
                                    blk{k} = blk_temp;
                                    
                                    BB2 = B2(:,blk_temp);
                                    CC2 = C2(:,blk_temp);
                                    AA2 = A2(blk_temp,:);
                                    
                                    sig2 = sum(Tab.Block_Size(next_block_ind(k)));
                                    [sizB1,sizB2] = size(BB2);
                                    BB = kron(eye(sig2), BB2);
                                    if include_sigma_energy
                                        CC1 = CC2(:);
                                        CC = repmat(CC1, 1, sizB2 * sig2);
                                    else
                                        CC = kron(eye(sig2), CC2);
                                    end
                                    AA = kron(transpose(AA2),eye(sizB1));
                                    AI{k} = AA*(CC.*BB);
                                end
                                
                                %disp(['N: ', num2str(n_part), ' S: ', num2str(spin), ', dag: ', num2str(dag), ...
                                %    ', c_spin: ', num2str(c_spin), ', mu: ', num2str(site_index_mu), ', kappa: ', num2str(site_index_kappa)])
                                %tttemp = cell2mat(AI);
                                %disp(tttemp + conj(tttemp(index1_trans,:)))
                                
                                
                                A3 = A3 + cell2mat(AI);
                                
                            end
                            % global strategy (not necessary if block divided
                            %
                            % [~,sig2] = deal(sum(Tab.Block_Size(next_block_ind))); %should be squared
                            % [sizB1,sizB2] = size(B2);
                            % BB = kron(eye(sig2), B2);
                            % CC1 = C2(:);
                            % CC = repmat(CC1, 1, sizB2 * sig2);
                            % AA = kron(transpose(A2),eye(sizB1));
                            % A3 = AA*(CC.*BB);
                            
                            % if isempty(spin),
                            %    index_cut = cell2mat(sel_block(Qtab.Index('N', n_part + dag) ).');
                            % else
                            %     index_cut = cell2mat(sel_block(Qtab.Index('NS', n_part + dag, spin + dag*c_spin) ).');
                            % end
                            
                            
                            
                            %K(index1,index3) = K(index1,index3) + A3 + conj(A3(index1_trans,:));
                            
                        end
                        
                        
                    end
                    if numel(energy_subindex) >= 1
                        index3 = Tab.Vector_Index(next_block_ind);
                        
                        
                        % %K(index1,index3) = K(index1,index3) + A3(:,index_cut) + conj(A3(index1_trans,index_cut));
                        % NOTE: TAKE CARE: the first index has to be the second argument in mesh grid!
                        
                        
                        
                        super_operator_A = A3 + conj(A3(:,:));
                        L_A = abs(super_operator_A(:)) > 0;
                        
                        if sum(L_A) > 0
                            
                            [i2, i1] = meshgrid(index3, index1);
                            i1 = i1(L_A); %NOTE: Problem with row instead of column can only happen if numel of one index is 1
                            i2 = i2(L_A);
                            w = super_operator_A(L_A);
                            %                         sparse_index = max(sparse_index) + (1:sum(L_A));
                            %                         W(sparse_index) = super_operator_A(L_A);
                            %                         I1(sparse_index) = i1(L_A);
                            %                         I2(sparse_index) = i2(L_A);
                            
                            [I1,I2,W, sparse_index] = add2sparse(I1,I2,W,sparse_index, i1(:), i2(:), w(:));
                        end
                        
                        %sparse_index = max(sparse_index) + (1:numel(i1));
                        %W(sparse_index) = super_operator_A(:);
                        
                        %I1(sparse_index) = i1(:);
                        %I2(sparse_index) = i2(:);
                        %TODO: figure out which is correct!!!
                        %super_operator_A = A3 + conj(A3(index1_trans,:));
                        
                    end
                    
                    
                end
                
            end
            %%
            
        end
        
        super_operator_B = - X1 - conj(X1(index1_trans, index1_trans));
        
        L_B = abs(super_operator_B(:)) > 0;
        
        
        if sum(L_B) > 0
            [i2, i1] = meshgrid(index1, index1);
            
            %         sparse_index = max(sparse_index) + (1:sum(L_B));
            %         W(sparse_index) = super_operator_B(L_B);
            %         I1(sparse_index) = i1(L_B);
            %         I2(sparse_index) = i2(L_B);
            
            [I1,I2,W, sparse_index] = add2sparse(I1,I2,W,sparse_index, i1(L_B),i2(L_B), super_operator_B(L_B));
            
        end
        
        %K(index1, index1) =  K(index1, index1) - X1 - conj(X1(index1_trans, index1_trans)); %?
        
        % TODO: optimization for secular approximation possible!!!
        %K(index1, index1) =  K(index1, index1) - kron(eye(Tab.Block_Size(n_block)), X1) - kron(conj(X1), eye(Tab.Block_Size(n_block))); %?
        if ~secular_approximation
            [Ea, Eb] = meshgrid(Er);
            
            [i2, i1] = meshgrid(index1, index1);
            super_operator_C =  - 1i*diag(Ea(:) - Eb(:));
            
            %         sparse_index = max(sparse_index) + (1:numel(i1));
            %         I1(sparse_index) = i1(:);
            %         I2(sparse_index) = i2(:);
            %         W(sparse_index) = super_operator_C(:);
            
            [I1,I2,W, sparse_index] = add2sparse(I1,I2,W,sparse_index, i1(:),i2(:),super_operator_C(:));
            %K(index1, index1) = K(index1, index1) - 1i*diag(Ea(:) - Eb(:));
        end
        sigma_index = [sigma_index; index1(:)];
        
        %TODO:  remove later (just for testing)
        
        
        Tab.Categ = Tab.Categ.f_recreate_lists;
        Tab_cut = Tab.f_cut_filter('I', energ_sort_index(1:n_block_index));
        Tab_cut.Energy = [];
        diag_index_shortened = Tab_cut.Diagonal_Index();
        [K, column_sum(n_block_index)]= sparse2K(I1,I2,W, sigma_index, max(sparse_index), diag_index_shortened, Tab_cut, mes_flag);
        
        
        
        ev = eig(full(K));
        ev_min(n_block_index) = min(abs(ev));
        % make video of convergence!
        %        figure, plot(real(ev), imag(ev), 'kx'), grid on, title(['block id: ', num2str(n_block_index)])
        [~, choi(n_block_index), choi_ew] = Tab_cut.f_choi(K);
       
        [Sig, status] = get_sigma(K, Tab_cut, status, 0, Energy_cut, Energy_cut_2, sigma_tolerance);
        if ~isempty(Sig)
            Sig_result{n_block_index} = Sig.s_dom_blocks;
        end
        
        
        if any(n_block_index == break_set)
            %check if K submits a stationary solution
            status.energy_cut = max(Er);
            col_sum_index = col_sum_index +1;
            Tab.Categ = Tab.Categ.f_recreate_lists;
            Tab_cut = Tab.f_cut_filter('I', energ_sort_index(1:n_block_index));
            Tab_cut.Energy = [];
            diag_index_shortened = Tab_cut.Diagonal_Index();
            
            %increase tolerance
            if mod(col_sum_index, 10) == 0
                status.eigenvalue_tolerance = status.eigenvalue_tolerance * 10;
            end
            
            % first and second:  generate K and reduce K to relevant matrix
            % third: calculate column sum
            [K, status.column_sum, column_sum_detail{col_sum_index}]= sparse2K(I1,I2,W, sigma_index, max(sparse_index), diag_index_shortened, Tab_cut, mes_flag);
            col_sum(col_sum_index) = status.column_sum;
            
            if size(K,1) < 5000
                ev = eig(full(K));
                figure, plot(real(ev), imag(ev), 'kx'), grid on, title(['block id: ', num2str(n_block_index)])
            else
                ev = eigs(K,3,10^-10);
            end
            ev_temp = ev(abs(ev) == min(abs(ev)));
            ev_vec(col_sum_index) = ev_temp(1);
            if status.column_sum < 0.1 %NOTE: add dynamic column_sum tolerance
                
                [Sig, status] = get_sigma(K, Tab_cut, status, mes_flag, Energy_cut, Energy_cut_2, sigma_tolerance);
                
                [~, choi_herm, choi_ew] = Tab_cut.f_choi(K);
                status.choi_herm = choi_herm;
                status.choi_neg_ew_ratio = sum(real(choi_ew) < -10^-10)/sum(~isnan(choi_ew));
                status.choi_mean_ew = mean(real(choi_ew(~isnan(choi_ew))));
                disp(['number of nan elements: ', num2str(sum(isnan(choi_ew))/numel(choi_ew))])
                if mes_flag
                    disp(status)
                end
                %status.choi =
                
                
                if ~status.no_solution %& status.choi_herm < 10^2
                    final_block_index = n_block_index;
                    disp(['-----------final eigenvalue tolerance: ', num2str(status.eigenvalue_tolerance)])
                    disp(['-----------minimal eigenvalue: ', num2str(status.minimal_value)])
                    disp(['-----------final_block_index: ', num2str(final_block_index), ' of ', num2str(numel(energ_sort_index))])
                    run_K_flag = false;
                else
                    if condor.save_intermediate_K
                        K_data = vars2struct('K', 'Tab_cut', 'Sig', 'status', 'Energy_cut', 'Energy_cut_2');
                        K_file_name = ['.', SLASH, 'log_files', SLASH, 'K_', condor.setup_id, '_', run_id, '_', num2str(n_block_index), '.mat'];
                        save(K_file_name, 'K_data')
                    end
                    
                    if final_block_index == n_block_index
                        if length(K) < 3000
                            ew = eig(full(K));
                            plot(real(ew), imag(ew), 'x')
                            grid on
                            title(['minimal value: ', num2str(min(abs(ew))), ', eigenvalue tolerance: ', num2str(status.eigenvalue_tolerance)])
                        end
                        error(['Convergence Error: no solution found - reset energy cut: ', num2str(Energy_cut)])
                        
                    end
                    
                    %continue
                end
            end
            
            % if column sum is sufficiently small, calculate eigenvalue zero
            % if found break, else continue to next break
            
        end
    end
    
end


%K = load('./log_files/K_0_1x1_430.mat')
%TODO: reprogram current
%%
disp('============== end of K generation ===================')
I = zeros(N_bath, En.Basis.N_spin);
%TODO: take only those blocks into account which are dominant (Sig.s_dom_blocks)
limit = 10^-10;
relevant_Table = Sig.s_dom_blocks(limit);
old_indices = relevant_Table.Old_Indices;
new_indices = relevant_Table.Index;

for n_block_index = 1:size(relevant_Table,1)
    
    rho = Sig.Values{new_indices(n_block_index)};
    
    for alpha = 1:N_bath
        
        for c_spin = c_spin_range
            
            I(alpha,c_spin/2 + 1.5) = I(alpha,c_spin/2 + 1.5) + 2*real(trace(X{old_indices(n_block_index), alpha, c_spin/2+1.5} * rho));
            
            
        end
    end
end


if abs(sum(I(:))) > 10^-6
    status.current_missmatch = sum(I(:));
end




