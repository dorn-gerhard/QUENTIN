function [final_nu_c, final_nu_m, final_nu_t, final_nu_trace,  final_energy, final_block_number, final_block_size, correlation, final_convergence_time, final_error, final_purity]  = ...
          check_convergence_2(K, Tab_cut, calc_time, Sig, status, Energy_cut, Energy_cut_2, start_block, En, tolerance)
    
if nargin < 10 || isempty(tolerance), tolerance = 10^-7; end
if nargin < 8 || isempty(start_block), start_block = ones(1, numel(tolerance)); end


sigma_tolerance = 10^-13;

if ismember('A', Tab_cut.Order)
    [energ_sort, energ_sort_index] = sort(Tab_cut.Average_Energies);
else
    [energ_sort, energ_sort_index] = sort(Tab_cut.Energies);
end



numel_tolerance = numel(tolerance);
index_tolerance = 1;
run_flag = true;

[final_nu_c, final_nu_m, final_nu_t, final_nu_trace, final_energy, final_block_number, final_block_size, correlation, final_convergence_time, final_error, final_purity] = deal(inf(1,numel_tolerance));

for block_index = 1:numel(energ_sort_index)
    if run_flag && block_index >= start_block(index_tolerance)
        
        % get reduced K_temp and Tab_temp
        vec_index = Tab_cut.Vector_Index('I', energ_sort_index(1:block_index));
        K_temp = K(vec_index, vec_index);
        
        % try to cut a cut table - problems with Li_Index and Energy!
        Tab_cut_temp = Tab_cut;
        Tab_cut_temp.Categ = Tab_cut_temp.Categ.f_remove('O');
        
        Tab_temp = Tab_cut_temp.f_cut_filter('I', energ_sort_index(1:block_index));
        
        
        % get minimal eigenvalue
        if size(K_temp,1) < 5000
            ev = eig(full(K_temp));
            %figure, plot(real(ev), imag(ev), 'kx'), grid on, title(['block id: ', num2str(n_block_index)])
        else
            ev = eigs(K_temp,3,10^-10);
        end
        ev_temp = ev(abs(ev) == min(abs(ev)));
        
        
        if abs(ev_temp) < tolerance(index_tolerance)
            
            final_nu_c(index_tolerance) = abs(ev_temp);
            [final_nu_m(index_tolerance), final_nu_t(index_tolerance)] = markovianity_trace(K_temp, Tab_temp);
            
            
            
            % get column sum
            K_relevant = K_temp(Tab_temp.Diagonal_Index(),:);
            %renormalize according to K-size
            final_nu_trace(index_tolerance) = sum(abs(sum(K_relevant,1)),2) / (size(K_relevant,2) * size(K_relevant,1));
            
            
            
            final_energy(index_tolerance) = energ_sort(block_index);
            
            
            final_block_number(index_tolerance) = block_index;
            
            
            final_block_size(index_tolerance) = sum(Tab_cut.Block_Size(energ_sort_index(1:block_index)).^2); % == length(K_temp)
            
            status.eigenvalue_tolerance = tolerance(index_tolerance);
            
            % get Sigma and compare with final sigma (vector norm?) diagonals more important? define a weight!!!
            [Sig_temp] = get_sigma(K_temp, En, Tab_temp, status, 0, Energy_cut, Energy_cut_2, sigma_tolerance);
                       
            
            Sig_mat = Sig_temp.f_matrix;
            Sig_mat_offdiag = Sig_mat - diag(diag(Sig_mat));
            correlation(index_tolerance) = (sum(abs(Sig_mat_offdiag(:))));
            
            final_purity(index_tolerance) = denoise(trace(Sig_mat^2), sigma_tolerance);
            
            final_convergence_time(index_tolerance) = calc_time(block_index);
            
            index_matrix_compare = Tab_cut.Matrix_Index(Tab_temp.Old_Indices);
            Sig_compare = sparse(index_matrix_compare(:,1), index_matrix_compare(:,2), Sig_temp.f_vector, Tab_cut.Numel_Datasets, Tab_cut.Numel_Datasets);
            DIFF = Sig_compare - Sig.f_matrix;
            
            final_error(index_tolerance) = sqrt(sum(abs(DIFF(:)).^2));
            
            
            
            
            if index_tolerance == numel_tolerance
                run_flag = false;
            else
                index_tolerance = index_tolerance + 1;
            end
            
            
            
        end
    end
    
end
