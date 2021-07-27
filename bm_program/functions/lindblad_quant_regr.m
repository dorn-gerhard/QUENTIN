function [K,I,Sig, Bat, Tab_cut, status, calc_time, column_sum, ev_min, choi, X_cur, H_LS_calculated] = lindblad_quant_regr(En, contact_site, Tab, Mu, Coupling, Beta, Bat, condor, Energy_cut, Energy_cut_2, run_id, min_en_ind, eigenvalue_tolerance, multi_symmetrization)

%%

if ismember('S', Tab.Order); spin_order = 'S'; else, spin_order = []; end
%Qtab = Table(En, ['N', spin_order]);
if nargin < 14 || isempty(multi_symmetrization), multi_symmetrization = false; end
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
% Qtab.Category(k).Value = Qtab.Category(k).Value(Lnew);
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


%NOTE: maybe preinitialize by max block number?
[column_sum, ev_min, choi, calc_time] = deal(0);


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
    energ_index = Tab.Block_Index(table_tag, [], [], {-inf, Energy_cut});
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
break_set = [min_en_ind:10:numel(energ_sort_index), final_block_index]; %final_block_index; %
break_set = 10:1:final_block_index;
col_sum = zeros(size(break_set));
col_sum_index = 0;

X_cur = cell( numel(energ_sort_index) , N_bath, numel(c_spin_range));
%X_sym = cell( numel(energ_sort_index) , N_bath, numel(c_spin_range));
%Bath3D = cell(N_bath,1);
%Bath3D_2 = cell(N_bath,1);
Bath3D_G = cell(N_bath,1);
Bath3D_S = cell(N_bath,1);

run_K_flag = true;
%energ_sort_index = energ_sort_index([1])
%%
H_LS_calculated = cell(final_block_index,1);

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
        
        [A2, H_LS] = deal(zeros(block_size));
        
        
        
        index1_trans = reshape(1:block_size^2, block_size,block_size)';
       
        
        % NOTE: more efficient (check if also fast)
        list_index = Tab.List_Index('I', n_block);
        
        %think about if Average Energy would be better?
        E_bd = En.Energies('I', list_index);
        sub_block_index = En.Subindex('NS',  {n_part, spin}, 'I', {list_index});
        numel_bd = numel(E_bd);      
        [X_cur{n_block, :, :}] = deal(0);

        for dag = -1:2:1
            
            for c_spin = c_spin_range %for spinless c_spin_range is just -1
                %%
                c_spin_index = c_spin/2+1.5;
                if spinless
                    c_spin = [];
                    c_spin_index = 1;
                end
                if numel(Tab.Block_Index(table_tag, n_part-dag, spin - dag*c_spin, {-inf, Energy_cut_2})) >= 1

                    count = count +1;
                    %NOTE: Aim of this section is to calculate Gamma_{ab,cd} and to calculate the entries in K
                    %Therefore get
                    %- Energies E_bd and E_ac
                    %- Q matrices
                    %- couplings
                    %- gamma and sigma terms (depending on energy differences)
                    
                    %energy_subindex = En.f_energy_subindex(n_part + dag, spin + dag*c_spin, {-inf, Energy_cut});
                    next_block_sizes = Tab.Block_Size(table_tag, n_part - dag, spin - dag*c_spin, {-inf, Energy_cut_2});
                    next_block_ind = Tab.Block_Index(table_tag, n_part - dag, spin - dag*c_spin, {-inf, Energy_cut_2});
                    
                    
                    
                    %getting the energies
                    E_ac = En.Energies(['NS', energy_tag], n_part - dag, spin - dag*c_spin, {-inf, Energy_cut_2});
                    [E_AC, E_BD] = meshgrid( E_ac, E_bd ); % creates matrices of the form < Enr | Enl> (N x N-s)
                    numel_ac = numel(E_ac);
                    
                    %getting Q matrices                    
                    energy_subindex = En.Subindex('NS', {n_part - dag, spin - dag*c_spin}, energy_tag, {{-inf, Energy_cut}});
                    energy_subindex_virtual = En.Subindex('NS',{n_part - dag, spin - dag*c_spin}, energy_tag, {{-inf, Energy_cut_2}});
                    %TODO: CHECK whether two subindexes (two energy cuts) are needed!
                    A_mu = zeros(numel_bd, numel_ac, numel(rel_site));

                    for site_index_mu = 1:numel(rel_site)
                        c_site_mu = rel_site(site_index_mu);
                        A_full = En.g_Qmat(dag, n_part-dag, spin-dag*c_spin, c_site_mu, c_spin, Energy_cut_2);
                        A_mu(:,:, site_index_mu) = full(A_full(sub_block_index, energy_subindex_virtual));
                    end
%                     for site_index_kappa = 1:numel(rel_site)
%                         c_site_kappa = rel_site(site_index_kappa);
%                         B_full = En.g_Qmat(-dag,n_part,spin, c_site_kappa, c_spin, Energy_cut_2);
%                         B_kappa(1:numel_ac, 1:numel_bd, site_index_kappa) = full(B_full(energy_subindex_virtual,sub_block_index));
%                     end
                    B_kappa = permute(conj(A_mu), [2,1,3]);
                    
                    if multi_symmetrization, mul_sym_fac = 1/2; else, mul_sym_fac = 1; end
                    %Bath3D_G = zeros(N_Orbitals, N_Orbitals, numel_bd, numel_ac, N_bath);
                    
                    
                    
                    %getting couplings
                    C_mu = permute(A_mu, [3,4,2,1]); 
                    C_kappa = permute(conj(A_mu), [4,3,5,6,1,2]);
                    
                    %getting couplings, gamma and sigma terms and Gamma (abcd)
                    
                    if numel(rel_site)^2 * numel_ac^2 * numel_bd^2 > 10^8
                        
                        % ================== work everything in blocks (next blocks) ===================
                        % ======================================================================================
                        % ======================================================================================
                        
                        
                        Gamma = cell(numel(next_block_ind) ,1);
                        A1_cell = cell(numel(next_block_ind) ,1);
                        for n_next_block = 1:numel(next_block_ind) 
                            next_block = next_block_ind(n_next_block);
                            sub_index = (1:next_block_sizes(n_next_block)) + sum(next_block_sizes(1:n_next_block-1));
                            numel_next_block = next_block_sizes(n_next_block);
                            
                            Gamma{n_next_block} = zeros(numel_next_block, numel_bd, numel_next_block, numel_bd);
                            
                            
                            for alpha = 1:N_bath
                                mu_alpha = Mu(alpha, 1);
                                beta_alpha = Beta(alpha,1);
                                V = Coupling{alpha,1}(rel_site,:);
                                %TODO: if V is nonzero for just one site, implement a saving strategy (1 instead of
                                %N_site^2 storage needed!
                                
                                Bath3D_G{alpha} = reshape(Bat.f_gamma(mul_sym_fac*(reshape(E_BD(:, sub_index)-E_AC(:, sub_index),numel_bd * numel_next_block,[] )), mu_alpha, beta_alpha, alpha, c_spin, dag), Bat.Numel_Orbitals(alpha), Bat.Numel_Orbitals(alpha), numel_bd, numel_next_block);
                                [Sigma_temp, Bat] =  Bat.f_lookup_sigma(mul_sym_fac*(reshape(E_BD(:, sub_index)-E_AC(:, sub_index),numel_bd * numel_next_block,[] )), mu_alpha, beta_alpha, alpha, c_spin, dag);
                                
                                %compare =  reshape(Bat.f_sigma(mul_sym_fac*(E_BD(:)-E_AC(:)), mu_alpha, beta_alpha, alpha, c_spin, dag), sqrt(Bat.Numel_Orbitals(alpha)),sqrt(Bat.Numel_Orbitals(alpha)),size(E_AC,1),size(E_AC,2));
                                
                                Bath3D_S{alpha} = reshape(Sigma_temp, Bat.Numel_Orbitals(alpha),Bat.Numel_Orbitals(alpha), numel_bd, numel_next_block);
                                
                                %Bath3D_S{alpha} = compare;
                                %NOTE: Bath3D_{alpha}(Orbital_l, Orbital_k, E_bd, E_ac) is Hermitian (conj(permute(,[2,1,...])
                                %NOTE: for working gmdmp: G(Orbital_l, E_bd, E_ac, Orbital_k) is AntiHermitian
                                G = permute(Bath3D_G{alpha}, [1,3,4,2]);
                                S = permute(Bath3D_S{alpha}, [1,3,4,2]);
                                
                                %Coupling{alpha, Bat.Spin_fun(c_spin)}' * Bath3D_G{alpha} * Coupling{alpha,Bat.Spin_fun(c_spin)}
                                
                                Gam = gmdmp(gmdmp(V, 2, 2, G, 1, 4), 4, 4, V', 1, 2); %sum over Orbitals k and l and baths \alpha
                                Sig = gmdmp(gmdmp(V, 2, 2, S, 1, 4), 4, 4, V', 1, 2);
                                
                                %--------------- calculate \Gamma_{ab,cd} ----------------
                                
                                % T (mu_index, kappa_index, a_index, b_index, d_index, c_index)
                                
                                
                                if sum(abs(Gam(:)) > 10^-12) > 0
                                    if numel_next_block > 300
                                        blk_size = 200;
                                        Gam_cur = zeros(numel_next_block, numel_bd, numel_next_block, numel_bd);
                                        blocks = mat2cell(1:numel_next_block, 1,[repmat(blk_size, 1, floor(numel_next_block/blk_size)), rem(numel_next_block, blk_size)]);
                                        for blocks_1_index = 1:numel(blocks)
                                            for blocks_2_index = 1:numel(blocks)
                                                blocks_1 = blocks{blocks_1_index};
                                                blocks_2 = blocks{blocks_2_index};
                                                if multi_symmetrization
                                                    T = 2*(permute(Gam(:,:,blocks_1,:), [1,4,3,2]) .* permute(Gam(:,:,blocks_2,:), [1,4,5,6,2,3]) + 0 * permute(Sig(:,:,blocks_1,:), [1,4,3,2]) .* permute(Sig(:,:,blocks_2,:), [1,4,5,6,2,3]));
                                                else
                                                    T = 1/2*(permute(Gam(:,:,blocks_1,:), [1,4,3,2]) + permute(Gam(:,:,blocks_2,:), [1,4,5,6,2,3]));
                                                end
                                                Gam_cur(blocks_1, :, blocks_2,:) = permute(sum(sum(C_mu(:,:,sub_index(blocks_1),:).*T.*C_kappa(:,:,:,:,:,sub_index(blocks_2)),1),2), [3,4,6,5,1,2]);
                                                clear T
                                            end
                                            
                                        end
                                    else
                                        if multi_symmetrization
                                            T = 2*(permute(Gam, [1,4,3,2]) .* permute(Gam, [1,4,5,6,2,3]) + 0 * permute(Sig, [1,4,3,2]) .* permute(Sig, [1,4,5,6,2,3]));
                                        else
                                            T = 1/2*(permute(Gam, [1,4,3,2]) + permute(Gam, [1,4,5,6,2,3]));
                                        end
                                        Gam_cur = permute(sum(sum(C_mu(:,:,sub_index,:).*T.*C_kappa(:,:,:,:,:,sub_index),1),2), [3,4,6,5,1,2]);
                                        clear T
                                    end
                                else
                                    Gam_cur = zeros(numel_next_block, numel_bd, numel_next_block, numel_bd);
                                    
                                    
                                end
                                
                                
                                if sum(abs(Sig(:)) > 10^-12) > 0
                                    if numel_next_block > 300
                                        blk_size = 200;
                                        Pi = zeros(numel_next_block, numel_bd, numel_next_block, numel_bd);
                                        blocks = mat2cell(1:numel_next_block, 1,[repmat(blk_size, 1, floor(numel_next_block/blk_size)), rem(numel_next_block, blk_size)]);
                                        for blocks_1_index = 1:numel(blocks)
                                            for blocks_2_index = 1:numel(blocks)
                                                blocks_1 = blocks{blocks_1_index};
                                                blocks_2 = blocks{blocks_2_index};
                                                if multi_symmetrization
                                                    T_LS = (permute(Gam(:,:,blocks_1,:), [1,4,3,2]) .* permute(Sig(:,:,blocks_2,:), [1,4,5,6,2,3]) +  permute(Sig(:,:,blocks_1,:), [1,4,3,2]) .* permute(Gam(:,:,blocks_2,:), [1,4,5,6,2,3])) ;
                                                else
                                                    T_LS = 1/4*(permute(Sig(:,:,blocks_1,:), [1,4,3,2]) + permute(Sig(:,:,blocks_2,:), [1,4,5,6,2,3]));
                                                end
                                                Pi(blocks_1,:,blocks_2,:) = permute(sum(sum(C_mu(:,:,sub_index(blocks_1),:).*T_LS.*C_kappa(:,:,:,:,:,sub_index(blocks_2)),1),2), [3,4,6,5,1,2]);
                                                clear T_LS
                                            end
                                            
                                        end
                                    else
                                        if multi_symmetrization
                                            T_LS = (permute(Gam, [1,4,3,2]) .* permute(Sig, [1,4,5,6,2,3]) +  permute(Sig, [1,4,3,2]) .* permute(Gam, [1,4,5,6,2,3])) ;
                                        else
                                            T_LS = 1/4*(permute(Sig, [1,4,3,2]) + permute(Sig, [1,4,5,6,2,3]));
                                        end
                                        Pi = permute(sum(sum(C_mu(:,:,sub_index,:).*T_LS.*C_kappa(:,:,:,:,:,sub_index),1),2), [3,4,6,5,1,2]);
                                        clear T_LS
                                    end
                                    
                                else 
                                    Pi = zeros(numel_next_block, numel_bd, numel_next_block, numel_bd);
                                    
                                end
                               
                                
                                X_cur{n_block, alpha, c_spin_index} = X_cur{n_block, alpha, c_spin_index} - dag * diagsum(Gam_cur,1,3);
                                Gamma{n_next_block} = Gamma{n_next_block} + Gam_cur;
                                
                                % Lamb Shift
                                
                                % Pi_abcd = \sum_{mu\kappa}  C_{ab,kappa}^{\overline{s}} \sum_\alpha V_kappa^s
                                % (g_{labk}^\alpha * s_{lcdk}^\alpha + s_{labk}^\alpha *g_{lcdk}^\alpha) * V_\mu^\overline{s}
                                % C_{dc\mu}^s

                                H_LS = H_LS  - 1i*permute(diagsum(Pi, 1,3),[2,1]);
                            end
                            
                            
                            %--------------------- construction of superoperator K from Gamma ----------------------------
                        
                            % \Gamma_{abcd} |a > < b| \sigma |d > < c|
                            A1_cell{n_next_block} = reshape(permute(Gamma{n_next_block}, [1,3,2,4]), numel_next_block^2, numel_bd^2);
                            
                            % \Gamma_{abcd} {|d > < c| |a > < b|, \sigma}
                            % A2 * \sigma + \sigma * A2
                            A2 = A2 - 1/2* permute(diagsum(Gamma{n_next_block}, 1,3),[2,1]);
                            
                            %TTest{n_next_block} = reshape(permute(Gamma{n_next_block}, [2,1,4,3]), numel_next_block * numel_bd, []);
                            
                            
                           
                            index2 = Tab.Vector_Index(next_block_ind(n_next_block));
                            
                            
                            L_A = abs(A1_cell{n_next_block}(:)) > 10^-14; %denoising
                            
                            if sum(L_A) > 0
                                
                                [i1, i2] = meshgrid(index1, index2);
                                i1 = i1(L_A); %NOTE: Problem with row instead of column can only happen if numel of one index is 1
                                i2 = i2(L_A);
                                w = A1_cell{n_next_block}(L_A);
                                
                                %NOTE: index1 is index of sigma on the right side, index2 on the left side (\dot{\sigma})
                                [I1,I2,W, sparse_index] = add2sparse(I1,I2,W,sparse_index, i2(:), i1(:), w(:));
                            end
  
                        end
                        
  
                    else
                        % ================== work everything in one step (small next blocks) ===================
                        % ======================================================================================
                        % ======================================================================================
                        
                        Gamma = zeros(numel_ac, numel_bd, numel_ac, numel_bd);
                        
                        for alpha = 1:N_bath
                            mu_alpha = Mu(alpha, 1);
                            beta_alpha = Beta(alpha, 1);
                            V = Coupling{alpha,1}(rel_site,:);
                            Bath3D_G{alpha} = reshape(Bat.f_gamma(mul_sym_fac*(E_BD(:)-E_AC(:)), mu_alpha, beta_alpha, alpha, c_spin, dag), Bat.Numel_Orbitals(alpha),Bat.Numel_Orbitals(alpha),size(E_AC,1),size(E_AC,2));
                            [Sigma_temp, Bat] =  Bat.f_lookup_sigma(mul_sym_fac*(E_BD(:)-E_AC(:)), mu_alpha, beta_alpha, alpha, c_spin, dag);
                            
                            %compare =  reshape(Bat.f_sigma(mul_sym_fac*(E_BD(:)-E_AC(:)), mu_alpha, beta_alpha, alpha, c_spin, dag), sqrt(Bat.Numel_Orbitals(alpha)),sqrt(Bat.Numel_Orbitals(alpha)),size(E_AC,1),size(E_AC,2));
                            
                            Bath3D_S{alpha} = reshape(Sigma_temp, Bat.Numel_Orbitals(alpha),Bat.Numel_Orbitals(alpha),size(E_AC,1),size(E_AC,2));
                            
                            %Bath3D_S{alpha} = compare;
                            %NOTE: Bath3D_{alpha}(Orbital_l, Orbital_k, E_bd, E_ac) is Hermitian (conj(permute(,[2,1,...])
                            %NOTE: for working gmdmp: G(Orbital_l, E_bd, E_ac, Orbital_k) is AntiHermitian
                            G = permute(Bath3D_G{alpha}, [1,3,4,2]);
                            S = permute(Bath3D_S{alpha}, [1,3,4,2]);
                            
                            %Coupling{alpha, Bat.Spin_fun(c_spin)}' * Bath3D_G{alpha} * Coupling{alpha,Bat.Spin_fun(c_spin)}
                            
                            Gam = gmdmp(gmdmp(V, 2, 2, G, 1, 4), 4, 4, V', 1, 2); %sum over Orbitals k and l and baths \alpha
                            Sig = gmdmp(gmdmp(V, 2, 2, S, 1, 4), 4, 4, V', 1, 2);
                            
                            %--------------- calculate \Gamma_{ab,cd} ----------------
                            
                            % T (mu_index, kappa_index, a_index, b_index, d_index, c_index)
                            
                            
                            
                            if multi_symmetrization %positive symmetric
                                T = 2* permute(Gam, [1,4,3,2]) .* permute(Gam, [1,4,5,6,2,3]) +0* permute(Sig, [1,4,3,2]) .* permute(Sig, [1,4,5,6,2,3]) ;
                                T_LS = (permute(Gam, [1,4,3,2]) .* permute(Sig, [1,4,5,6,2,3]) +  permute(Sig, [1,4,3,2]) .* permute(Gam, [1,4,5,6,2,3])) ;
                            else %symmetric
                                T =  1/2*(permute(Gam, [1,4,3,2]) + permute(Gam, [1,4,5,6,2,3]));
                                T_LS = 1/4*(permute(Sig, [1,4,3,2]) + permute(Sig, [1,4,5,6,2,3]));
                            end
                            
                            Gam_cur = permute(sum(sum(C_mu.*T.*C_kappa,1),2), [3,4,6,5,1,2]);%sum over kappa and mu
                            %Gamma(a,b,c,d) (switch from abdc to abcd)
                            clear T
                            
                            X_cur{n_block, alpha, c_spin_index} = X_cur{n_block, alpha,c_spin_index} - dag * diagsum(Gam_cur,1,3);
                            
                            Gamma = Gamma + Gam_cur;
                            
                            % H_LS (Lamb shift)
                            Pi = permute(sum(sum(C_mu.*T_LS.*C_kappa,1),2), [3,4,6,5,1,2]);
                            clear T_LS
                            H_LS = H_LS  - 1i*permute(diagsum(Pi, 1,3),[2,1]);
                            % final: Gam(mu_index, E_bd, E_ac, kappa_index)
                        end
                        
                        
                        %Test Gamma
                        %                     Super = reshape(Gamma, [numel_ac* numel_bd, numel_ac*numel_bd]);
                        %                     Herm_super = Super - Super';
                        %                     ew_super = eig(1/2*(Super+Super'));
                        %                     disp(['hermiticity of Gamma: ', num2str(sum(abs(Herm_super(:)))), ', number of negative eigenvalues: ', ...
                        %                         num2str(sum(real(ew_super) < -10^-10))])
                        %                     ew_super(real(ew_super) < -10^-10);
                        %                     ew_super(real(ew_super) > 10^-10);
                        %
                        %--------------------- construction of superoperator K from Gamma ----------------------------
                        
                        % \Gamma_{abcd} |a > < b| \sigma |d > < c|
                        %SS = reshape(permute(Gamma, [2,1,4,3]), numel_bd*numel_ac, numel_bd*numel_ac);
                        A1 = reshape(permute(Gamma, [1,3,2,4]), numel_ac^2, numel_bd^2);
                        
                        
                        % \Gamma_{abcd} {|d > < c| |a > < b|, \sigma}
                        
                        % A2 * \sigma + \sigma * A2
                        A2 = A2 - 1/2* permute(diagsum(Gamma, 1,3),[2,1]);
                        

                        % ------------------------------inserting in K
                        
                        
                        %NOTE: in order to cut quasi energy blocks from A1! (alternative: loop over blocks! ->
                        %implemented in born_markov_full_dynamic.m)
                        block_list = cell2mat(arrayfun(@(x,y) repmat(x,y,1), next_block_ind, next_block_sizes, 'uniformoutput', false));
                        [xx,yy] = meshgrid(block_list);
                        L_block = xx(:) == yy(:);
                        A1 = A1(L_block, :);
                        
                        index2 = Tab.Vector_Index(next_block_ind);
                        
                        
                        L_A = abs(A1(:)) > 10^-14; %denoising
                        
                        if sum(L_A) > 0
                            
                            [i1, i2] = meshgrid(index1, index2);
                            i1 = i1(L_A); %NOTE: Problem with row instead of column can only happen if numel of one index is 1
                            i2 = i2(L_A);
                            w = A1(L_A);
                            
                            %NOTE: index1 is index of sigma on the right side, index2 on the left side (\dot{\sigma})
                            [I1,I2,W, sparse_index] = add2sparse(I1,I2,W,sparse_index, i2(:), i1(:), w(:));
                        end
                        
                        
                    end
                    
                    % ======================================================================================
                    % ======================================================================================
                    % ======================================================================================
                    
                end
                
            end
            %%
            
        end
        
        
        
        
        super_operator_B = (kron(eye(block_size), A2) + kron(transpose(A2), eye(block_size)));
        
        L_B = abs(super_operator_B(:)) > 10^-14; %denoising
        
        
        if sum(L_B) > 0
            [i2, i1] = meshgrid(index1, index1);
            
            [I1,I2,W, sparse_index] = add2sparse(I1,I2,W,sparse_index, i1(L_B),i2(L_B), super_operator_B(L_B));
            
        end
        H_LS_calculated{n_block} = H_LS ;
        H_LS = H_LS + diag(E_bd);
        %Temp = -1i*(kron(eye(block_size), H_LS) - kron(transpose(H_LS), eye(block_size)));
        
        
        super_operator_C = -1i*(kron(eye(block_size), H_LS) - kron(transpose(H_LS), eye(block_size)));
        
        L_C = abs(super_operator_C(:)) > 10^-14; %denoising
        
        
         if sum(L_C) > 0
            [i2, i1] = meshgrid(index1, index1);
            
            [I1,I2,W, sparse_index] = add2sparse(I1,I2,W,sparse_index, i1(L_C),i2(L_C), super_operator_C(L_C));
        end
        
        
        
        %K(index1, index1) = K(index1, index1) - 1i*diag(Ea(:) - Eb(:));
        
        
        
        
        sigma_index = [sigma_index; index1(:)];
        
        %TODO:  remove later (just for testing)
        
        if false
            check_choi = false;
            Tab.Categ = Tab.Categ.f_recreate_lists;
            Tab_cut = Tab.f_cut_filter('I', energ_sort_index(1:n_block_index));
            Tab_cut.Energy = [];
            diag_index_shortened = Tab_cut.Diagonal_Index();
            %NOTE(Index I1 and I2 flipped (cutoff in first index!)
            cut_off_first_index = true;
            [K, column_sum(n_block_index)]= sparse2K(I1,I2,W, sigma_index, max(sparse_index), diag_index_shortened, Tab_cut, mes_flag, cut_off_first_index);
            
            % make video of convergence!
            ev = eig(full(K));
            %figure, plot(real(ev), imag(ev), 'kx'), grid on, title(['block id: ', num2str(n_block_index)])
            if check_choi
                [~, choi(n_block_index), choi_ew] = Tab_cut.f_choi(K);
            else
                choi(n_block_index) = nan;
            end
            
            
            if mes_flag
                [Sig, status] = get_sigma(K, Tab_cut, status, 0, Energy_cut, Energy_cut_2, sigma_tolerance);
                ev_min(n_block_index) = status.minimal_value;
                disp(['column_sum: ', num2str(column_sum(n_block_index), 5), ', choi: ', num2str(choi(n_block_index),5), ', min_ev: ', num2str(ev_min(n_block_index),5)])
                
                if ~isempty(Sig)
                    Sig_result{n_block_index} = Sig.s_dom_blocks;
                    Sig.s_dom_blocks
                end
                
            end
        end
        
        calc_time(n_block_index) = toc(calc_time_counter);
        if any(n_block_index == break_set)
            %check if K submits a stationary solution
            status.energy_cut = max(E_bd);
            col_sum_index = col_sum_index +1;
            
            Tab.Categ = Tab.Categ.f_recreate_lists;
            Tab_cut = Tab.f_cut_filter('I', energ_sort_index(1:n_block_index)); %maintains the ordering of the original table
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
            cut_off_first_index = true;
            [K, status.column_sum, column_sum_detail{col_sum_index}]= sparse2K(I1,I2,W, sigma_index, max(sparse_index), diag_index_shortened, Tab_cut, mes_flag, cut_off_first_index);
            col_sum(col_sum_index) = status.column_sum;
            
%             if size(K,1) < 5000
%                 ev = eig(full(K));
%                 %figure, plot(real(ev), imag(ev), 'kx'), grid on, title(['block id: ', num2str(n_block_index)])
%             else
%                 ev = eigs(K,3,10^-10);
%             end
%             if ~isempty(ev)
%                 ev_temp = ev(abs(ev) == min(abs(ev)));
%                 ev_vec(col_sum_index) = ev_temp(1);
%             else
%                 ev_vec(col_sum_index) = inf;
%             end
            if status.column_sum < 0.1 %NOTE: add dynamic column_sum tolerance
                
                [Sig, status] = get_sigma(K, Tab_cut, status, mes_flag, Energy_cut, Energy_cut_2, sigma_tolerance);
                if ~status.no_solution
                    tab_full =  Sig.Tab.s_table('O', energ_sort_index(1:n_block_index));
                    if ismember('A', Tab_cut.Order)
                        
                        [energ_sort, Index_sort] = sort(Tab_cut.Average_Energies);
                    else
                        [energ_sort, Index_sort] = sort(Tab_cut.Energies);
                    end
                    weight = tab_full.Weight(Index_sort);
                    
                    %energ_sort > max(energ_sort(weight < status.eigenvalue_tolerance))
                    %if any(full(column_sum_detail{col_sum_index}(:)) < status.eigenvalue_tolerance & weight < status.eigenvalue_tolerance)
                    if column_sum_detail{col_sum_index}(Index_sort(end)) < status.eigenvalue_tolerance || any(full(column_sum_detail{col_sum_index}(:)) < status.eigenvalue_tolerance & energ_sort > max(energ_sort(weight >= status.eigenvalue_tolerance)))
                        %converged
                    else
                        status.no_solution = 1;
                        status.error_text = ['trace number of non_relevant states does not match tolerance:', num2str(status.eigenvalue_tolerance), ', minimal trace number: ' ...
                            , num2str(min(column_sum_detail{col_sum_index}(weight < status.eigenvalue_tolerance)))];
                    end
                end
              
                
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
            
            %NOTE: maybe reduce by factor of 2 since F = 1/2 ( gamma + sigma)
            
            I(alpha,c_spin/2 + 1.5) = I(alpha,c_spin/2 + 1.5) + real(sum(sum(X_cur{old_indices(n_block_index), alpha, c_spin/2+1.5} .* rho)));
            
            
        end
    end
end


if abs(sum(I(:))) > 10^-6
    status.current_missmatch = sum(I(:));
end




