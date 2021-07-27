function [C, length] = print_aligned_table(A, delimiter)
if nargin < 2 || isempty(delimiter), delimiter = '  '; end


%A = [12 0; 2.12 0.000212];
[siz_x,siz_y] = size(A);


tt = num2str(A(:), '%3.3f');

L = true(size(tt,1),1);
for k = 1:size(tt,2)
    if all(tt(L,end-k+1) == '.')
        L =  mod(A(:),1) == 0;
        
        tt(L,end-k+1) = ' ';
        break
    end
    L2 = tt(L,end-k+1) == '0';
    L(L) = L2;
    tt(L,end-k+1) = ' ';
    
    
end


tt = [tt, repmat(delimiter, size(tt,1),1)];
length = size(tt,2);
C = cell(size(tt,1),1);
for k = 1:size(tt,1)
    C{k,1} = tt(k,:);
end

C = reshape(C, siz_x,siz_y);


