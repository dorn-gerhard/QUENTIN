function [K,I,Sig, Bat, Tab_cut, status, calc_time, column_sum, ev_min, choi_herm, X_cur] = born_markov_full_dynamic(En, contact_site, Tab, Mu, Coupling, Beta, Bat, condor, Energy_cut, Energy_cut_2, run_id, min_en_ind, eigenvalue_tolerance, symmetrized_version)

%%
% sort the blocks by their energy
if ismember('A', Tab.Order)
    [energ_sort, energ_sort_index] = sort(Tab.Average_Energies);
else
    [energ_sort, energ_sort_index] = sort(Tab.Energies);
end

if ismember('S', Tab.Order); spin_order = 'S'; else, spin_order = []; end
%Qtab = Table(En, ['N', spin_order]);
if nargin < 14 || isempty(symmetrized_version), symmetrized_version = true; end
if nargin < 13 || isempty(eigenvalue_tolerance), eigenvalue_tolerance = 10^-7; end
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
status.partial_trace_converged = inf;
status.energy_domain_limit = inf;
status.steady_state_negative = inf;
status.choi_herm = inf;

status.nu_m = inf;
status.nu_t = inf;
status.full_gamma_positivity = inf;

status.error_text = '';

[column_sum, ev_min, choi_herm, calc_time] = deal(0);


%SWITCH for Redfield Bloch type K_0
include_sigma_energy = false;


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



%[energ_sort, energ_sort_index]
energ_sort_index(energ_sort  > Energy_cut) = []; %delete blocks that exceed Energy_cut

final_block_index = numel(energ_sort_index);
break_set = [min_en_ind:10:final_block_index, final_block_index]; %final_block_index; %

col_sum = zeros(size(break_set));
col_sum_index = 0;

X_cur = cell( numel(energ_sort_index) , N_bath, numel(c_spin_range));
%X_sym = cell( numel(energ_sort_index) , N_bath, numel(c_spin_range));
%Bath3D = cell(N_bath,1);
%Bath3D_2 = cell(N_bath,1);
Bath3D_G1 = cell(N_bath,1);
Bath3D_G2 = cell(N_bath,1);
Bath3D_S = cell(N_bath,1);
Bath3D_S2= cell(N_bath,1);

run_K_flag = true;
%energ_sort_index = energ_sort_index([1])
%%

calc_time_counter = tic;
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
        XX = zeros(block_size);
        Y1 = zeros(block_size);
        
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
        E_N = En.Energies('I', list_index);
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
        [X_cur{n_block, :, :}] = deal(0);
%         [X_sym{n_block, :, :}] = deal(0);

        
        for dag = -1:2:1
            
            for c_spin = c_spin_range %for spinless c_spin_range is just -1
                %%
               c_spin_index = c_spin/2+1.5;
                if spinless
                    c_spin = [];
                    c_spin_index = 1;
                end
                
                count = count +1;
                
                %disp(['round: ', num2str(count)])
                %NOTE: querry for virtual excitations: a^s a^_s \rho  (Energy_cut_2) case A)
                %NOTE: for second part only Energy_cut relevant a^s \rho a^_s        case B)
                if numel(Tab.Block_Index(table_tag, n_part-dag, spin - dag*c_spin, {-inf, Energy_cut_2})) >= 1
                    
                    %TODO: Test with 'A' instead of E
                    energy_subindex = En.Subindex('NS', {n_part - dag, spin - dag*c_spin}, energy_tag, {{-inf, Energy_cut}});
                    
                    %energy_subindex = En.f_energy_subindex(n_part + dag, spin + dag*c_spin, {-inf, Energy_cut});
                    next_block_sizes = Tab.Block_Size(table_tag, n_part - dag, spin - dag*c_spin, {-inf, Energy_cut});
                    next_block_ind = Tab.Block_Index(table_tag, n_part - dag, spin - dag*c_spin, {-inf, Energy_cut});
                    
                    index2_trans = cell2mat(arrayfun(@(x,y) reshape(reshape(x+(1:y^2), y,[]).',[],1), cumsum(next_block_sizes.^2)-next_block_sizes.^2, next_block_sizes, 'uniformOutput', false));
                    
                    %TODO: Test with 'A' instead of E
                    energy_subindex_virtual = En.Subindex('NS',{n_part - dag, spin - dag*c_spin}, energy_tag, {{-inf, Energy_cut_2}});
                    %energy_subindex_virtual = En.f_energy_subindex(n_part + dag, spin + dag*c_spin, {-inf, Energy_cut_2});
                    
                    %TODO: Test with 'A' instead of E
                    
                    %There are two operators needed: A * B * sigma
                    
                    %getting the energies
                    E_M = En.Energies(['NS', energy_tag], n_part - dag, spin - dag*c_spin, {-inf, Energy_cut_2});
                    [E_MM, E_NN] = meshgrid( E_M, E_N ); % creates matrices of the form < Enr | Enl> (N x N-s)
                    
                    E_M2 = En.Energies(['NS', energy_tag], n_part - dag, spin - dag*c_spin, {-inf, Energy_cut_2});
                    [E_MM2, ~] = meshgrid( E_M2, E_N ); % creates matrices of the form < Enr | Enl> (N x N-s)
                    
                    A3 = zeros(block_size^2 , sum(next_block_sizes.^2) ); %relevant for case B)
                    
                    
                    [G1,G2,S,S2] = deal(zeros(numel(E_N), numel(E_M), En.Basis.N_site, En.Basis.N_site, N_bath));
                    for alpha = 1:N_bath
                        mu_alpha = Mu(alpha, 1);
                        beta_alpha = Beta(alpha, 1);
                        
                        if dag == 1
                            coupling_mu = conj(Coupling{alpha,1});
                            coupling_kappa = Coupling{alpha,1};
                        else
                            coupling_mu = Coupling{alpha,1};
                            coupling_kappa = conj(Coupling{alpha,1});
                        end
                        % first term, form: A*(G1.*B)*sigma; A = c_\mu^s, B = c_\kappa^{\overline{s}}
                        
                        
                            
                        
                        Bath3D_G1{alpha} = reshape(Bat.f_gamma((E_NN(:)-E_MM(:)), mu_alpha, beta_alpha, alpha, c_spin, dag), (Bat.Numel_Orbitals(alpha)),(Bat.Numel_Orbitals(alpha)),size(E_MM,1),size(E_MM,2));
                        G1(:,:,:,:,alpha) = permute(mmat(mmat(coupling_mu,Bath3D_G1{alpha}), coupling_kappa'),[3,4,1,2]); %indexing: Er, El, mu, kappa, alpha
                        
                        %G2 = \gamma_{\kappa\mu}^{\overline{s}} (E_(N) - E_N-s), the change of s has to be
                        %considered in the couplings, in the dagger; NOTE that the energy has a different sign,
                        %also exchange \kappa and \mu
                        Bath3D_G2{alpha} = reshape(Bat.f_gamma(-(E_NN(:) - E_MM2(:)), mu_alpha, beta_alpha, alpha, c_spin, -dag), (Bat.Numel_Orbitals(alpha)),(Bat.Numel_Orbitals(alpha)),size(E_MM,1),size(E_MM,2));
                        G2(:,:,:,:,alpha) = permute(mmat(mmat(conj(coupling_kappa),Bath3D_G2{alpha}), conj(coupling_mu')),[3,4,1,2]); %indexing: Er, El, mu, kappa, alpha
                        
                        
                        [Sigma_temp, Bat] = Bat.f_lookup_sigma((E_NN(:)-E_MM(:)), mu_alpha, beta_alpha, alpha, c_spin, dag);
                        Bath3D_S{alpha} =  reshape(Sigma_temp, (Bat.Numel_Orbitals(alpha)),(Bat.Numel_Orbitals(alpha)),size(E_MM,1),size(E_MM,2));
                        S(:,:,:,:,alpha) = permute(mmat(mmat(coupling_mu,Bath3D_S{alpha}), coupling_kappa'),[3,4,1,2]); %indexing: Er, El, mu, kappa, alpha
                        
                        if ~symmetrized_version
                            [Sigma_temp_2, Bat] = Bat.f_lookup_sigma(-(E_NN(:)-E_MM2(:)), mu_alpha, beta_alpha, alpha, c_spin, -dag);

                            Bath3D_S2{alpha} =  reshape(Sigma_temp_2, (Bat.Numel_Orbitals(alpha)),(Bat.Numel_Orbitals(alpha)),size(E_MM,1),size(E_MM,2));
                            S2(:,:,:,:,alpha) = permute(mmat(mmat(conj(coupling_kappa),Bath3D_S2{alpha}), conj(coupling_mu')),[3,4,1,2]); %indexing: Er, El, mu, kappa, alpha
                        end
                        
                        
%                         Bath3D_S{alpha} =  reshape(Bat.f_sigma(mul_sym_fac*(E_NN(:)-E_MM(:)), mu_alpha, beta_alpha, alpha, c_spin, dag), sqrt(Bat.Numel_Orbitals(alpha)),sqrt(Bat.Numel_Orbitals(alpha)),size(E_MM,1),size(E_MM,2));
%                         S(:,:,:,:,alpha) = permute(mmat(mmat(coupling_mu,Bath3D_S{alpha}), coupling_kappa'),[3,4,1,2]); %indexing: Er, El, mu, kappa, alpha
%                         
%                         if ~symmetrized_version
%                             Bath3D_S2{alpha} =  reshape(Bat.f_sigma(-mul_sym_fac*(E_NN(:)-E_MM2(:)), mu_alpha, beta_alpha, alpha, c_spin, -dag), sqrt(Bat.Numel_Orbitals(alpha)),sqrt(Bat.Numel_Orbitals(alpha)),size(E_MM,1),size(E_MM,2));
%                             S2(:,:,:,:,alpha) = permute(mmat(mmat(conj(coupling_kappa),Bath3D_S2{alpha}), conj(coupling_mu')),[3,4,1,2]); %indexing: Er, El, mu, kappa, alpha
%                         end
%                         
%                         
                        
                        
                        
                        
                    end
                    
                    %size(G2)
                    
                    % ================ check Gammas of Lindblad form ============
                    if false
                        % create Gamma matrix in that energy block
                        % blockwise!!
                        for k = 1:numel(next_block_ind)
                            
                            blk_temp = (1:next_block_sizes(k)) + sum(next_block_sizes(1:k))-next_block_sizes(k);
                            G2_temp = sum(G2(:,blk_temp,:,:,:),5);
                            %orders first the energies, then the sites (E1,site1|E1, site1), (E1,site2| E1,site1), ...
                            TT = permute(G2_temp,[3,1,4,2]);
                            [numel_EN, numel_EM, size_mu, size_kappa] = size(G2_temp);
                            TT3 = reshape(TT,size_mu * numel_EN,size_kappa * numel_EM);
                            
                            TT4 = reshape(permute(G2_temp,[3,2,4,1]), size_mu * numel_EM,size_kappa * numel_EN);
                            
                            %GAM is hermitian!
                            %check smallest eigenvalue
                            if size(TT3,1) == size(TT3,2)
                                GAM =  TT3 + TT4;
                                GAM_ew = real(eig(GAM));
                            else
                                GAM_ew = real(eig(sqrt(TT3*TT4)));
                            end
                            if min(GAM_ew) < -10^-5
                                disp(['dag: ', num2str(dag), ', c_spin: ', num2str(c_spin), ', minew: ', num2str(min(GAM_ew)), ', max: ', num2str(max(GAM_ew))])
                            end
                            
                            %TT2 = [squeeze(G2_temp(1,1,:,:)), squeeze(G2_temp(1,2,:,:)); squeeze(G2_temp(2,1,:,:)), squeeze(G2_temp(2,2,:,:))]
                            
                            %ew = eig([squeeze(G2_temp(1,1,:,:)), squeeze(G2_temp(1,2,:,:)); squeeze(G2_temp(2,1,:,:)), squeeze(G2_temp(2,2,:,:))])
                        end
                    end
                    
                    for site_index_mu = 1:numel(rel_site)
                        c_site_mu = rel_site(site_index_mu);
                        
                        
                        
                        A_full = En.g_Qmat(dag, n_part-dag, spin-dag*c_spin, c_site_mu, c_spin, Energy_cut_2);
                        A = A_full(sub_block_index, energy_subindex_virtual);
                        
                        for site_index_kappa = 1:numel(rel_site)
                            c_site_kappa = rel_site(site_index_kappa);
                          
                            
                            
                            B_full = En.g_Qmat(-dag,n_part,spin, c_site_kappa, c_spin, Energy_cut_2);
                            B = B_full(energy_subindex_virtual,sub_block_index); % here opposite, since A is kind of adjoint
                            
                            %for the current calculation:
                            for alpha = 1:N_bath
                                if symmetrized_version %TODO: think about symmetrized current formula
                                    %X_sym{n_block, alpha, Bat.Spin(fun(c_spin))} = X_sym{n_block, alpha, Bat.Spin_fun(c_spin)} -1i/4*(A_virt*(gamma_B(:,:,c_site_mu,c_site_kappa, alpha).*B_virt) + (gamma_A(:,:,c_site_mu,c_site_kappa,alpha).*A_virt)*B_virt);

                                    X_cur{n_block, alpha, c_spin_index} =X_cur{n_block, alpha, c_spin_index} + dag * A*(G1(:,:, c_site_mu, c_site_kappa,alpha).'.*B);
                                else
                                    F2 = S(:,:, c_site_mu, c_site_kappa,alpha) + G1(:,:, c_site_mu, c_site_kappa,alpha);
                                    X_cur{n_block, alpha, c_spin_index} =X_cur{n_block, alpha, c_spin_index} + dag * A*( F2.'.*B);
                                end
                            end
                
                                                        
                            % ------------------first term---------------------
                            % form (A*(G1.'.*B) + (A.*G1)*B ) * sigma
                            % c_mu^s c_kappa^{\overline{s}} -> 
                            % A = < E  | c_mu^s | E'> , 
                            % B = < E' |c_kappa^{\overline{s}} | E >
                            % 
                            sig_size = Tab.Block_Size(n_block);
                            
                            GG1 = sum(G1(:,:,c_site_mu, c_site_kappa,:),5); %sum over baths
                            SS = sum(S(:,:,c_site_mu, c_site_kappa,:),5);
                            
                            if symmetrized_version
                                
                                X = A*(B.*GG1.') + (A.*GG1)*B;
                                Y = (SS.*A)*B + A*(SS.'.*B);
                                %disp(-1/4*X)
                                %disp(eig(-1/4*full(X)))
                                XX = XX + X;
                                
                                
                                %Lamb shift term:
                                
                                Y1 = Y1 + Y;
                            else
                                X = 2* A*(B.*(GG1.' + SS.'));
                               
                            end
                            %NOTE: the 1/4 stems from F = 1/2 (gamma + sigma) in both cases
                            X1 = X1 - 1/4 * kron(eye(sig_size),X);
                            
                            % % for inclusion of sigma
                            % BB = kron(eye(sig_size), B_virt);
                            % CC_B = kron(eye(sig_size), gamma_B(:,:,c_site_mu, c_site_kappa,alpha));
                            %
                            % AA = kron(eye(sig_size), A_virt);
                            % CC_A = kron(eye(sig_size), gamma_A(:,:,c_site_mu, c_site_kappa,alpha));
                            %
                            % if symmetrized_version
                            %   X1 = X1 + 1/2 * (AA*(CC_B.*BB) + (AA.*CC_A) * BB);
                            % else
                            % 	X1 = X1 + AA * (CC_B .* BB);
                            % end
                            if false
                                [~,sig_size] = deal(Tab.Block_Size(n_block)); %should be squared
                                
                                [~,sizB2] = size(B_virt);
                                C = sum(C,3);
                                C_sym = sum(C_sym,3);
                               
                                BB = kron(eye(sig_size), B_virt);
                                
                                
                                
                                if include_sigma_energy
                                    CC1 = C(:);
                                    CC = repmat(CC1, 1, sizB2 * sig_size);
                                else
                                    CC = kron(eye(sig_size), C);
                                    CC_sym = kron(eye(sig_size),C_sym);
                                end
                                AA = kron(eye(sig_size), A_virt);
                                
                                %disp([num2str(n_part),', S: ', num2str(spin),', dag: ',  num2str(dag), ', c_spin: ' , num2str(c_spin)])
                                %disp(['N: ', num2str(n_part), ' S: ', num2str(spin), ', dag: ', num2str(dag), ...
                                %    ', c_spin: ', num2str(c_spin), ', mu: ', num2str(site_index_mu), ', kappa: ', num2str(site_index_kappa)])
                                
                                %disp(AA*(CC.*BB))
                                if symmetrized_version
                                    X1 = X1 + 1/2*( AA*(CC.*BB) + (AA.*CC_sym)*BB);
                                else
                                    X1 = X1 + AA*(CC.*BB); % diagonal block element
                                end
                            end
                            %AA*(CC.*BB)
                            % TODO: think of complex conjugate
                            
                            
                            
                            % ------------------second term-------------------
                            
                            % form (G2.*A) * sigma * B + h.c.
                            % (\gamma_{\kappa\mu}^{overline{s}} c_\mu^s) \sigma c_\kappa^{\overline{s}} ->
                            % A = < E  | c_mu^s | E'> ,
                            % B = < E' |c_kappa^{\overline{s}} | E >
                            % 
                            
                            if numel(energy_subindex) >= 1
                                %TODO: reduce energies which are in between Energy_cut and Energy_cut_2
                            
                                
                                A2 = A_full(sub_block_index, energy_subindex);
                                B2 = B_full(energy_subindex,sub_block_index);
                                %cut also G2: 
                                if symmetrized_version
                                    GG2 = sum(G2(:,:,c_site_mu, c_site_kappa,:),5);
                                else
                                    GG2 = sum(G2(:,:,c_site_mu, c_site_kappa,:),5) + sum(S2(:,:,c_site_mu, c_site_kappa,:),5);
                                end
                                GG2_tilde = GG2(:,energy_subindex);
                                                                
                                
                                
                                %next_block_ind lists the blocks in the next particle / spin sector with energies
                                %below energy_cut
                                AI = cell(1,numel(next_block_ind));
                              
                                for k = 1:numel(next_block_ind)
                                    blk_temp = (1:next_block_sizes(k)) + sum(next_block_sizes(1:k))-next_block_sizes(k);

                                    BB2 = B2(blk_temp,:);
                                   
                                    
                                    
                                    sig_size = Tab.Block_Size(next_block_ind(k));
                                    [sizB2,sizB1] = size(BB2);
                                    
                                    BB = kron(transpose(BB2),eye(sizB1));
                                    
                                    if include_sigma_energy
                                        AA = kron(eye(sig_size), B2(:,blk_temp)); %B without
                                        CC1 = reshape(GG2(:,energy_subindex),[],1);
                                        CC = repmat(CC1, 1, sizB2 * sig_size);
                                        AI{k} = BB*(CC.*AA);
                                    else
                                        %NOTE: Test to symmetrize (make it like a diagonal) in order to get
                                        %positivity, maybe there is a better way.
                                        
                                        GG2_sym = GG2_tilde(:,blk_temp);
                                        
                                        
                                        AA2 = GG2_sym .* A2(:,blk_temp);
                                        
                                        
                                        AA = kron(eye(sig_size), AA2);
                                        AI{k} = BB*AA;
                                    end
                                    
                                    
                                end
                                
                                %disp(['N: ', num2str(n_part), ' S: ', num2str(spin), ', dag: ', num2str(dag), ...
                                %    ', c_spin: ', num2str(c_spin), ', mu: ', num2str(site_index_mu), ', kappa: ', num2str(site_index_kappa)])
                                %tttemp = cell2mat(AI);
                                %disp(tttemp + conj(tttemp(index1_trans,:)))
                                
                                %The 1/2 stems from the F = (gamma + sigma)/2
                                A3 = A3 + 1/2*cell2mat(AI); % A3 + 1/2 * cell2mat(AI);
                                
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
                        index2 = Tab.Vector_Index(next_block_ind);
                        
                        
                        % %K(index1,index3) = K(index1,index3) + A3(:,index_cut) + conj(A3(index1_trans,index_cut));
                        % NOTE: TAKE CARE: the first index has to be the second argument in mesh grid!
                       
                        
                        
                        super_operator_A = A3 + conj(A3(index1_trans,index2_trans));
                        L_A = abs(super_operator_A(:)) > 10^-14; %denoising
                        
                        if sum(L_A) > 0
                            
                            [i2, i1] = meshgrid(index2, index1);
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
        
        super_operator_B = X1 + conj(X1(index1_trans, index1_trans));
        
        L_B = abs(super_operator_B(:)) > 10^-14; %denoising
        
        
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
        
        [Ea, Eb] = meshgrid(E_N);
        [i2, i1] = meshgrid(index1, index1);
        super_operator_C =  - 1i*diag(Eb(:) - Ea(:)); %(left minus right)
        if symmetrized_version
            super_operator_C = super_operator_C - 1/4*(kron(eye(block_size), Y1) - kron(transpose(Y1), eye(block_size)));
        end
        
        %         sparse_index = max(sparse_index) + (1:numel(i1));
        %         I1(sparse_index) = i1(:);
        %         I2(sparse_index) = i2(:);
        %         W(sparse_index) = super_operator_C(:);
        
        [I1,I2,W, sparse_index] = add2sparse(I1,I2,W,sparse_index, i1(:),i2(:),super_operator_C(:));
        %K(index1, index1) = K(index1, index1) - 1i*diag(Ea(:) - Eb(:));
        
        
        
        
        sigma_index = [sigma_index; index1(:)];
        
        %TODO:  remove later (just for testing)
        
        if false
            Tab.Categ = Tab.Categ.f_recreate_lists;
            Tab_cut = Tab.f_cut_filter('I', energ_sort_index(1:n_block_index));
            Tab_cut.Energy = [];
            diag_index_shortened = Tab_cut.Diagonal_Index();
            [K, column_sum(n_block_index)]= sparse2K(I1,I2,W, sigma_index, max(sparse_index), diag_index_shortened, Tab_cut, mes_flag);
            
            % make video of convergence!
            ev = eig(full(K));
            % figure, plot(real(ev), imag(ev), 'kx'), grid on, title(['block id: ', num2str(n_block_index)])
            if check_choi
                [~, choi(n_block_index), choi_ew] = Tab_cut.f_choi(K);
            else
                choi(n_block_index) = nan;
            end
            
            
            if mes_flag
                [Sig, status] = get_sigma(K, En, Tab_cut, status, 0, Energy_cut, Energy_cut_2, sigma_tolerance);
                ev_min(n_block_index) = status.minimal_value;
                disp(['column_sum: ', num2str(column_sum(n_block_index), 5), ', choi: ', num2str(choi(n_block_index),5), ', min_ev: ', num2str(ev_min(n_block_index),5)])
                
               
                
            end
        end
        
        calc_time(n_block_index) = toc(calc_time_counter);
        if any(n_block_index == break_set)
            %check if K submits a stationary solution
            status.energy_cut = max(E_N);
            col_sum_index = col_sum_index +1;
            
            Tab.Categ = Tab.Categ.f_recreate_lists;
            Tab_cut = Tab.f_cut_filter('I', energ_sort_index(1:n_block_index));
            Tab_cut.Energy = [];
            diag_index_shortened = Tab_cut.Diagonal_Index();
            
            %increase tolerance
            if ~isfield(condor, 'increase_tolerance')
                condor.increase_tolerance = 10;
            end
            if mod(col_sum_index, condor.increase_tolerance) == 0
                status.eigenvalue_tolerance = status.eigenvalue_tolerance * 10;
            end
            
            % first and second:  generate K and reduce K to relevant matrix
            % third: calculate column sum
            cut_off_first_index = false;
            [K, status.column_sum, column_sum_detail{col_sum_index}]= sparse2K(I1,I2,W, sigma_index, max(sparse_index), diag_index_shortened, Tab_cut, mes_flag, cut_off_first_index);
            %NOTE: column_sum_detail is sorted by Average_Energies!!!
            col_sum(col_sum_index) = status.column_sum;
            
%             if size(K,1) < 5000
%                 ev = eig(full(K));
%                 %figure, plot(real(ev), imag(ev), 'kx'), grid on, title(['block id: ', num2str(n_block_index)])
%             else
%                 ev = eigs(K,10,10^-10);
%             end
%             ev_temp = ev(abs(ev) == min(abs(ev)));
%             ev_vec(col_sum_index) = ev_temp(1);
            if status.column_sum < 0.1 %NOTE: add dynamic column_sum tolerance
                
                [status.nu_m, status.nu_t] = markovianity_trace(K, Tab_cut);
                
                [Sig, status] = get_sigma(K, En, Tab_cut, status, mes_flag, Energy_cut, Energy_cut_2, sigma_tolerance);
                
                
                
                choi_herm = Tab_cut.f_choi_analyse(K);
                status.choi_herm = choi_herm;
                %status.choi_neg_ew_ratio = sum(real(choi_ew) < -10^-10)/sum(~isnan(choi_ew));
                %status.choi_mean_ew = mean(real(choi_ew(~isnan(choi_ew))));
                %disp(['number of nan elements: ', num2str(sum(isnan(choi_ew))/numel(choi_ew))])
                if mes_flag
                    disp(status)
                end
%                 %status.choi =
%                 
                
                if ~status.no_solution %& status.choi_herm < 10^2
                    status.steady_state_negative = Sig.Analysis.sum_neg_ev;
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
            
            %NOTE: maybe reduce by factor of 2 since F = 1/2 ( gamma + sigma)
            
            I(alpha,c_spin/2 + 1.5) = I(alpha,c_spin/2 + 1.5) - real(trace(X_cur{old_indices(n_block_index), alpha, c_spin/2+1.5} * rho));
            
            
        end
    end
end


if abs(sum(I(:))) > 10^-6
    status.current_missmatch = sum(I(:));
end




