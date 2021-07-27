function [all_cycle, description] = f_all_cycle(A)

n = size(A,1);

B = cell(0);
for k = 1:n
    B = [B, f_cycle(A,k, [], k)];
end


min_ind = cellfun(@(x) find(x == min(x)), B);
for k = 1:length(B)
    if min_ind(k) ~= 1
        B{k} = [B{k}(min_ind(k):end), B{k}(1:min_ind(k)-1)];
    end
end

L_one_direction = cellfun(@(x) x(2) < x(end), B);

description = unique(cellfun(@mat2str, B(L_one_direction), 'uniformoutput', false))';


all_cycle = cellfun(@str2num, description, 'uniformoutput', false);
       
    