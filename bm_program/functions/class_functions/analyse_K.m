function [blocks, ev_vec, choi_herm, column_sum, column_sum_detail, block_size, energ_sort, sig_comp_norm, weight_comp] = ...
    analyse_K(K, Tab_cut, calc_time, Sig, status, Energy_cut, Energy_cut_2, ID, K_ID, print_now)



sigma_tolerance = 10^-10;

if ismember('A', Tab_cut.Order)
    [energ_sort, energ_sort_index] = sort(Tab_cut.Average_Energies);
else
    [energ_sort, energ_sort_index] = sort(Tab_cut.Energies);
end

table_final = Sig.s_dom_blocks;

[ev_vec, pos_eigenvalues, column_sum, choi_herm, sig_comp_norm, weight_comp] = deal(zeros(numel(energ_sort_index,1)));




col_sum_index = 0;
for k = 1:numel(energ_sort_index)
    col_sum_index = col_sum_index + 1;
    
    % get reduced K_temp and Tab_temp
    vec_index = Tab_cut.Vector_Index('I', energ_sort_index(1:k));
    K_temp = K(vec_index, vec_index);
    if false
        Tab.Categ = Tab.Categ.f_recreate_lists;
        Tab_temp = Tab.f_cut_filter('I', Tab_cut.Old_Indices(energ_sort_index(1:k)));
    else
        % try to cut a cut table - problems with Li_Index and Energy!
        Tab_cut_temp = Tab_cut;
        Tab_cut_temp.Categ = Tab_cut_temp.Categ.f_remove('O');
        
        Tab_temp = Tab_cut_temp.f_cut_filter('I', energ_sort_index(1:k));
    end
    
    % get minimal eigenvalue
    if size(K_temp,1) < 5000
        ev = eig(full(K_temp));
        %figure, plot(real(ev), imag(ev), 'kx'), grid on, title(['block id: ', num2str(n_block_index)])
    else
        ev = eigs(K_temp,3,10^-10);
    end
    ev_temp = ev(abs(ev) == min(abs(ev)));
    ev_vec(col_sum_index) = ev_temp(1);
    
    % get number of positive eigenvalues
    pos_eigenvalues(col_sum_index) = sum(real(ev)> 10^-10);
    
    
    
    % get column sum
    K_relevant = K_temp(Tab_temp.Diagonal_Index(),:);
    %renormalize according to K-size
    column_sum(col_sum_index) = sum(abs(sum(K_relevant,1)),2) / (size(K_relevant,2) * size(K_relevant,1));
  
    
    if ismember('A', Tab_cut.Order)
        
        [Ener_sort, Index_sort] = sort(Tab_temp.Average_Energies);
    else
        [Ener_sort, Index_sort] = sort(Tab_temp.Energies);
    end

    
    for r = 1:k
        column_sum_detail{col_sum_index}(r) = sum(abs(sum(K_relevant(:,Tab_temp.Vector_Index(Index_sort(1:r))),1)),2) / (numel(Tab_temp.Vector_Index(Index_sort(1:r)))* size(K_relevant,1));
    end
    
    
    
    % get choi number
    [~, choi_herm(col_sum_index), ~] = Tab_temp.f_choi(K_temp);
    %[~, choi_herm(col_sum_index), choi_ew] = Tab2.f_choi(K_temp);
    
    tolerance = 10^-6;
    status.eigenvalue_tolerance = tolerance;
    
    
    % get Sigma and compare with final sigma (vector norm?) diagonals more important? define a weight!!!
    [Sig_temp{col_sum_index}, status_temp{col_sum_index}] = get_sigma(K_temp, Tab_temp, status, 0, Energy_cut, Energy_cut_2, sigma_tolerance);
    %[Sig_temp, status_temp] = get_sigma(K_temp, Tab2, status, mes_flag, Energy_cut, Energy_cut_2, sigma_tolerance);
    
    if ~isempty(Sig_temp{col_sum_index})
        
       trans_index = Tab_cut.List_Index(Sig_temp{col_sum_index}.Tab.Old_Indices) ;
       Sig_trans = zeros(size(Sig.f_matrix));
       Sig_trans(trans_index, trans_index) = Sig_temp{col_sum_index}.f_matrix;
       
       Sig_comp = Sig.f_matrix - Sig_trans;
       
       sig_comp_norm(col_sum_index) = sqrt(sum(abs(Sig_comp(:)).^2));
       
       table_temp = Sig_temp{col_sum_index}.s_dom_blocks;
       [~, old_index] = ismember(table_temp.Old_Indices, table_final.Index);
       L = old_index ~= 0;
       if sum(old_index ~=0) > 0
           weight_comp(col_sum_index) = sum(abs(table_temp.weight(L) - table_final.weight(old_index(L))));
       else
           weight_comp(col_sum_index) = nan;
       end
       

    end
    
    disp(['block nr: ', num2str(k), '/', num2str(numel(energ_sort_index)), ', column_sum: ', num2str(column_sum(col_sum_index), 5), ', choi: ', num2str(choi_herm(col_sum_index),5), ', min_ev: ', num2str(ev_vec(col_sum_index),5)])

end
blocks = 1:numel(energ_sort_index);
block_size = Tab_cut.Block_Size(energ_sort_index); 

ev_vec = [ev_vec, zeros(1, numel(energ_sort_index) - numel(ev_vec))];
choi_herm = [choi_herm, zeros(1, numel(energ_sort_index) - numel(choi_herm))];
column_sum = [column_sum, zeros(1, numel(energ_sort_index) - numel(column_sum))];



fig = figure;
tabgrp = uitabgroup(fig);

conv_K_tab = uitab(tabgrp, 'Title', 'K convergence', 'backgroundcolor', 'w');
tabgrp.SelectedTab = conv_K_tab;
axes('parent', conv_K_tab)


subplot(3,1,1)
semilogy(calc_time, abs(ev_vec), '.-')
grid on
hold on
semilogy(calc_time, choi_herm, '--')
semilogy(calc_time, column_sum)
set(gcf,'color', 'w')
xlabel('calculation time / s', 'interpreter', 'Latex')
legend('steady state eigenvalue', 'choi number', 'column sum', 'interpreter', 'Latex')
title('Qualification of the superoperator generate with the ERMEA approach', 'interpreter', 'Latex')
ylim([10^-17, 1])

subplot(3,1,2)
yyaxis left

semilogy(blocks, abs(ev_vec), '.-')
grid on
hold on
semilogy(blocks, choi_herm, '--')
semilogy(blocks, column_sum)
xlabel('included block', 'interpreter', 'Latex')
ylim([10^-17, 1])

title('Qualification of the superoperator generate with the ERMEA approach', 'interpreter', 'Latex')
%TODO: add blocksize
yyaxis right
bar(blocks, Tab_cut.Block_Size(energ_sort_index), 'edgecolor', 'none')
alpha 0.5
ylim([0,max(3, max(Tab_cut.Block_Size))])
legend('steady state eigenvalue', 'choi number', 'column sum', 'block size', 'interpreter', 'Latex')

subplot(3,1,3)
semilogy(energ_sort, abs(ev_vec), '.-')
grid on
hold on
semilogy(energ_sort, choi_herm, '--')
semilogy(energ_sort, column_sum)
xlabel('included energy', 'interpreter', 'Latex')
title('Qualification of the superoperator generate with the ERMEA approach', 'interpreter', 'Latex')
plot_settings
colors =  colormap('lines');
for l = 1:5:numel(energ_sort_index)
    semilogy(energ_sort(1:l), column_sum_detail{l}, ':', 'color', colors(3,:))

end
legend('steady state eigenvalue', 'choi number', 'column sum', 'partial column sums', 'interpreter', 'Latex')

ylim([10^-17, 1])

plot_supertitle(['$\verb|', ID, ', K: ', num2str(K_ID), '|$'], conv_K_tab)

if print_now
    plot_pdf([ID, '_K_', num2str(K_ID)])
end

%plot_pdf(['convergence.pdf'])
% analyse Sigma convergence

conv_sig_tab = uitab(tabgrp, 'Title', '$\sigma$ convergence', 'backgroundcolor', 'w');
tabgrp.SelectedTab = conv_sig_tab;
axes('parent', conv_sig_tab)

subplot(3,1,1)
semilogy(calc_time, abs(sig_comp_norm), '.-')
grid on
hold on
set(gcf,'color', 'w')
xlabel('calculation time / s')
legend('Frobenius norm $\sigma_{\textnormal{end}} - \sigma_t$', 'convergence of weights', 'interpreter', 'Latex')
title('Comparison of intermediate reduced density matrix with final one')
semilogy(calc_time, abs(weight_comp), '--')
xlabel('calculation time / s', 'interpreter', 'Latex')
ylim([10^-17, 1])


subplot(3,1,2)
semilogy(blocks, abs(sig_comp_norm), '.-')
yyaxis left
hold on
grid on
semilogy(blocks, abs(weight_comp), '--')
xlabel('included block', 'interpreter', 'Latex')
ylim([10^-17, 1])

yyaxis right
bar(blocks, Tab_cut.Block_Size(energ_sort_index), 'edgecolor', 'none')
ylim([0,max(3, max(Tab_cut.Block_Size))])
alpha 0.5
legend('Frobenius norm $\sigma_{\textnormal{end}} - \sigma_t$', 'convergence of weights', 'block size', 'interpreter', 'Latex')



subplot(3,1,3)
semilogy(energ_sort, abs(sig_comp_norm), '.-')
hold on
grid on
semilogy(energ_sort, abs(weight_comp), '--')
xlabel('included energy', 'interpreter', 'Latex')
ylim([10^-17, 1])

plot_settings
plot_supertitle(['$\verb|', ID, ', K: ', num2str(K_ID), '|$'], conv_sig_tab)

if print_now
    plot_pdf([ID, '_K_', num2str(K_ID)])
end

end
