
function [K, column_sum, column_sum_detail] = sparse2K(I1,I2,W, sigma_index, max_sparse_index, diag_index_shortened, Tab, mes_flag, cut_off_first_index)
if nargin < 9 || isempty(cut_off_first_index), cut_off_first_index = false; end
I1 = I1(1:max_sparse_index,1);
I2 = I2(1:max_sparse_index,1);
W = W(1:max_sparse_index,1);
%K = sparse(I1,I2,W, sum(Tab.Block_Size.^2)  , sum(Tab.Block_Size.^2)  );

sort_sigma_index = sort(sigma_index);

if cut_off_first_index
    [~,II2] = ismember(I2,sort_sigma_index);
    [check_2, II1] =ismember(I1,sort_sigma_index);
else
    [~,II1] = ismember(I1,sort_sigma_index);
    [check_2, II2] =ismember(I2,sort_sigma_index);
end
% Cutoff in second index!!!

K = sparse(II1(check_2), II2(check_2), W(check_2), numel(sigma_index), numel(sigma_index));
K(abs(K)< 10^-15) = 0;


if  length(K) < 100000
     if ismember('E', Tab.Order)
         Ener = Tab.Energies;
     elseif ismember('A', Tab.Order)
         Ener = Tab.Average_Energies;
     end
    
    %      diag_index = Tab.Diagonal_Index('I', energ_index);
    %      [~,diag_index_shortened] = ismember(diag_index, sigma_index);
    K_relevant = K(diag_index_shortened,:);
    %renormalize according to K-size
    column_sum = sum(abs(sum(K_relevant,1)),2) / (size(K,2) * numel(diag_index_shortened));
  
    

    [Ener_sort, Index_sort] = sort(Ener);
    for k = 1:numel(Ener_sort)
        column_sum_detail(k) = sum(abs(sum(K_relevant(:,Tab.Vector_Index(Index_sort(1:k))),1)),2) / (numel(Tab.Vector_Index(Index_sort(1:k)))*numel(diag_index_shortened));
    end
    
    
end



end