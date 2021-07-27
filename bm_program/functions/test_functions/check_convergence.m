function [final_ev, final_choi, final_column_sum, column_sum_energy, final_energy, final_block_number, final_block_size, correlation, final_convergence_time, final_error, final_nu_m, final_nu_t] = ...
          check_convergence(K, Tab_cut, calc_time, Sig, status, Energy_cut, Energy_cut_2, start_block, En, tolerance)
if nargin < 9 || isempty(tolerance), tolerance = 10^-7; end



sigma_tolerance = 10^-10;

if ismember('A', Tab_cut.Order)
    [energ_sort, energ_sort_index] = sort(Tab_cut.Average_Energies);
else
    [energ_sort, energ_sort_index] = sort(Tab_cut.Energies);
end

table_final = Sig.s_dom_blocks;

[ev_vec, pos_eigenvalues, column_sum, choi_herm, sig_comp_norm, weight_comp, nu_m, nu_t] = deal(zeros(numel(energ_sort_index),1));
column_sum_detail = cell(numel(energ_sort_index),1);

converged = false;
final_column_sum = inf;
set_break = false;

for block_index = start_block:numel(energ_sort_index)
    
    
    % get reduced K_temp and Tab_temp
    vec_index = Tab_cut.Vector_Index('I', energ_sort_index(1:block_index));
    K_temp = K(vec_index, vec_index);
    if false
        Tab.Categ = Tab.Categ.f_recreate_lists;
        Tab_temp = Tab.f_cut_filter('I', Tab_cut.Old_Indices(energ_sort_index(1:block_index)));
    else
        % try to cut a cut table - problems with Li_Index and Energy!
        Tab_cut_temp = Tab_cut;
        Tab_cut_temp.Categ = Tab_cut_temp.Categ.f_remove('O');
        
        Tab_temp = Tab_cut_temp.f_cut_filter('I', energ_sort_index(1:block_index));
    end
    
    % get minimal eigenvalue
    if size(K_temp,1) < 5000
        ev = eig(full(K_temp));
        %figure, plot(real(ev), imag(ev), 'kx'), grid on, title(['block id: ', num2str(n_block_index)])
    else
        ev = eigs(K_temp,3,10^-10);
    end
    ev_temp = ev(abs(ev) == min(abs(ev)));
    ev_vec(block_index) = ev_temp(1);
    
    % get number of positive eigenvalues
    pos_eigenvalues(block_index) = sum(real(ev)> 10^-10);
    
    
    
    % get column sum
    K_relevant = K_temp(Tab_temp.Diagonal_Index(),:);
    %renormalize according to K-size
    column_sum(block_index) = sum(abs(sum(K_relevant,1)),2) / (size(K_relevant,2) * size(K_relevant,1));
  
    
    if ismember('A', Tab_cut.Order)
        
        [Ener_sort, Index_sort] = sort(Tab_temp.Average_Energies);
    else
        [Ener_sort, Index_sort] = sort(Tab_temp.Energies);
    end

    
    
    
    
    
    % get choi number
    choi_herm(block_index) = Tab_temp.f_choi_analyse(K_temp);
    % to check positivity with Choi matrix calculate the Choi matrix of expm(K) and check for positivity
    [nu_m(block_index), nu_t(block_index)] = markovianity_trace(K_temp, Tab_temp);
    
    
    
    if false
        [super_index, K2] = Tab_temp.f_choi_index_permutation(expm(K_temp));
        %[~, choi_herm(col_sum_index), choi_ew] = Tab2.f_choi(K_temp);
        
        ew = eig((full(1000*K_temp)))
        figure
        plot(real(ew), imag(ew), 'x')
        ew = eig(expm(full(K_temp)))
        hold on
        plot(real(ew), imag(ew), 'ro')
        
        for k = 1:1000
            t = rand()*10;
            S = Sigma([],[],Tab_temp);
            
            trace(Tab_temp.f_matrix(full(K_temp) * ((rand(10,1)-0.5) + 1i * (rand(10,1) - 0.5))))
            [sort(eig(full(S.f_matrix))), ...
            sort(real(eig(full(Tab_temp.f_matrix( expm(full(K_temp) * 700) * S.f_vector)))))]
            if all(eig(full(Tab_temp.f_matrix( expm(full(K_temp) * t) * S.f_vector))) > 0)
                %okay
                %disp('okay')
            else
                warning('negative eigenvalues')
            end
        end
        
        
    end
    status.eigenvalue_tolerance = tolerance;
    
    % get Sigma and compare with final sigma (vector norm?) diagonals more important? define a weight!!!
    [Sig_temp{block_index}, status_temp{block_index}] = get_sigma(K_temp, En, Tab_temp, status, 0, Energy_cut, Energy_cut_2, sigma_tolerance);
    %[Sig_temp, status_temp] = get_sigma(K_temp, Tab2, status, mes_flag, Energy_cut, Energy_cut_2, sigma_tolerance);
    
    if ~isempty(Sig_temp{block_index})
       
       trans_index = Tab_cut.List_Index(Sig_temp{block_index}.Tab.Old_Indices) ;
       Sig_trans = zeros(size(Sig.f_matrix));
       Sig_trans(trans_index, trans_index) = Sig_temp{block_index}.f_matrix;
       
       Sig_comp = Sig.f_matrix - Sig_trans;
       
       sig_comp_norm(block_index) = sqrt(sum(abs(Sig_comp(:)).^2));
       
       table_temp = Sig_temp{block_index}.s_dom_blocks;
       [~, old_index] = ismember(table_temp.Old_Indices, table_final.Index);
       if sum(old_index ~=0) > 0
           weight_comp(block_index) = sum(abs(table_temp.Weight - table_final.Weight(old_index)));
       else
           weight_comp(block_index) = nan;
       end
       

    end
    
    
    
    % get earliest convergence!!!
    
    if ~isempty(Sig_temp{block_index})
        tab_temp_full = Sig_temp{block_index}.Tab.s_table;
        for r = 1:block_index
            column_sum_detail{block_index}(r) = full(sum(abs(sum(K_relevant(:,Tab_temp.Vector_Index(Index_sort(1:r))),1)),2) / (numel(Tab_temp.Vector_Index(Index_sort(1:r)))* size(K_relevant,1)));
        end 
        if ismember('A', Sig_temp{block_index}.Tab.Order)
            [energ_sort, energ_sort_index_temp] = sort(Sig_temp{block_index}.Tab.Average_Energies);
        else
            [energ_sort, energ_sort_index_temp] = sort(Sig_temp{block_index}.Tab.Energies);
        end
        %NOTE: ignore choi herm at the moment, since for large matrices calcualtion takes a long time and it is
        %anyway always fullfilled
        L = find(full(column_sum_detail{block_index}(:)) < tolerance & energ_sort > max(energ_sort(tab_temp_full.Weight(energ_sort_index_temp) >= status.eigenvalue_tolerance)), 1, 'first');
        if ~isempty(L)
            final_column_sum = column_sum_detail{block_index}(L);
        end
        
        if abs(ev_vec(block_index)) < tolerance && any(full(column_sum_detail{block_index}(:)) < tolerance & energ_sort > max(energ_sort(tab_temp_full.Weight(energ_sort_index_temp) >= status.eigenvalue_tolerance)))
        
        %if abs(ev_vec(block_index)) < tolerance && any((column_sum_detail{block_index}(:) < tolerance) & (tab_temp_full.Weight(energ_sort_index_temp) < tolerance))
            
            if ismember('A', Sig_temp{block_index}.Tab.Order)
                column_sum_energy = tab_temp_full.Average_Energies(energ_sort_index_temp(L));
            else
                column_sum_energy = tab_temp_full.Energies(energ_sort_index_temp(L));
            end
            converged = true;
            % convergence reached return block number or other relevant parameters such as computation time - give
            % correspondence with final sigma!
            set_break = true;
        end
    end
    
    disp(['block nr: ', num2str(block_index), '/', num2str(numel(energ_sort_index)), ', column_sum: ', num2str(column_sum(block_index), 5), ', final_column_sum: ', num2str(final_column_sum), ', choi: ', num2str(choi_herm(block_index),5), ', min_ev: ', num2str(ev_vec(block_index),5), ...
        ', \sigma error: ', num2str(weight_comp(block_index))])
    
    if set_break
        break
    end
    
end
if ~converged
    [column_sum_detail{block_index}(:) < tolerance , tab_temp_full.Weight(energ_sort_index_temp) < tolerance]
    disp(['min ev: ', num2str(abs(ev_vec(block_index))), ', column_sum: ', num2str(min(column_sum_detail{block_index}(:)))])
end
blocks = 1:block_index;

ev_vec = [ev_vec];
choi_herm = [choi_herm];
column_sum = [column_sum];

final_block_size = length(K_temp);
Sig_mat = Sig_temp{block_index}.f_matrix;
Sig_mat = Sig_mat - diag(diag(Sig_mat));
correlation = sum(abs(Sig_mat(:)));
final_convergence_time = calc_time(block_index);
final_ev = ev_vec(block_index);
final_choi = choi_herm(block_index);
final_block_number = block_index;
final_energy = energ_sort(block_index);
final_error = weight_comp(block_index);
final_nu_m = nu_m(block_index);
final_nu_t = nu_t(block_index);


