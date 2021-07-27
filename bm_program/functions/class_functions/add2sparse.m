function [I1,I2,W, sparse_index] = add2sparse(I1,I2,W,sparse_index, i1,i2,w)

% I1 = [1,1,1,2,2,3,3,3,4,4,4,4,4]'; I2 = [1,2,3,2,4,1,4,5,1,2,4,5,6]';
% W = (1:13)'; sparse_index = 1:13; i1 = [2,2,3,3,3]'; i2 = [1,4,1,3,4]'; w = [1,2,3,4,5]';

if true
    
    I1_temp = I1(1:max(sparse_index));
    I2_temp = I2(1:max(sparse_index));
    W_temp = W(1:max(sparse_index));
    
    
%     [ii,~, icA] = unique([i1,i2], 'rows', 'sorted');
%     
%     [II, ~, IXC] = unique([I1_temp, I2_temp], 'rows', 'sorted');
%     
    
    rows = [I1_temp,I2_temp;i1,i2];
    index = sortrowsc(rows,[1,2]);
    
    d = rows(index(1:end-1),:) == rows(index(2:end),:); %double entries
    d = all(d,2);
    ndx1 = index(d); %indicates which indices are double in I1
    
    %indicates which indices are double in i1
    ndx2 = index([false;d]) - max(sparse_index);
    
    if numel(ndx1) > 0
        W_temp(ndx1) = W_temp(ndx1) + w(ndx2); %perform summation
    end
    
    
    
    lia = false(numel(i1), 1);
    lia(ndx2) = true;
    if numel(ndx2) < length(i1)
        % add new elements
        sparse_index = max(sparse_index) + (1:sum(~lia));
        I1_new = [I1_temp;i1(~lia)];
        I2_new = [I2_temp; i2(~lia)];
        W_new = [W_temp; w(~lia)];
        
        %resort
        index_new = sortrowsc([I1_new,I2_new], [1,2]);
        
        
        I1(1:max(sparse_index)) = I1_new(index_new);
        I2(1:max(sparse_index)) = I2_new(index_new);
        W(1:max(sparse_index)) = W_new(index_new);
    else
        
        W(1:max(sparse_index)) = W_temp;
        
        
    end
    
    
    
    


else
    sparse_index = max(sparse_index) + (1:numel(i1));
    I1(sparse_index) = i1;
    I2(sparse_index) = i2;
    W(sparse_index) = w;
end
%profile viewer
