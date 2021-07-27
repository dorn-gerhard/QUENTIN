function M = index2sparse(I1,I2,W, sigma_index, max_sparse_index, cut_off_first_index)
if nargin < 6 || isempty(cut_off_first_index), cut_off_first_index = false; end




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
W = denoise(W, 10^-18);
M = sparse(II1(check_2&II1), II2(check_2&II1), W(check_2&II1), numel(sigma_index), numel(sigma_index));


end

