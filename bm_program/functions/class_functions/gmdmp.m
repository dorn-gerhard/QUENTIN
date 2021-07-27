
function [gmdmproduct] = gmdmp(tensor1, d1, ndim1, tensor2, d2, ndim2)
%GMDMP General Multi Dimensional Matrix Product.
%C = GMDMP(A, d1, ndim1, B, d2, ndim2) Computes the product
%C(i[1],...,i[d1-1],i[d1+1],...,i[m],j[1],...,j[d2-1],j[d2+1],...,j[n]) =
%     A(i[1],...,i[d1-1],k,i[d1+1],...,i[m]) * B(j[1],...,j[d2-1],k,j[d2+1],...,j[n])
%(Sum on k).
%
%C = GMDMP(A, d1, B, d2) takes the outer product of A and B, then traces
%along the diagonal formed by dimensions d1 of A and d2 of B. For example,
%C = GMDMP(A, ndims(A), B, 1) is just the natural extension of 2D matrix
%multiplication, and for A and B both 2D, coincides with C = A * B.
%
%Note: it matters not if the lengths of dimensions d1 of A and d2 of B do
%not agree (see header of DIAGSUM.M).
%
%Wynton Moore, January 2006
%evaluate
ndim1_ = ndims(tensor1); 
ndim2_ = ndims(tensor2);
dims1 = [size(tensor1), ones(1, ndim1-ndim1_)]; 
dims2 = [size(tensor2), ones(1, ndim2-ndim2_)]; 
temp = outer(tensor1,tensor2,0);

if ndim1_ < ndim1
    resort = [1:ndim1_, ndim1_+ndim2 + (1:ndim1-ndim1_), ndim1_ + (1:ndim2)];
    temp = permute(temp, resort);
end

tot_dims = [dims1,dims2];
if tot_dims(d1) == 1 && tot_dims(d2 + ndim1) == 1
    perm = 1:numel(tot_dims); perm([d1,d2+ndim1]) = []; perm = [perm, d1,d2+ndim1];
    gmdmproduct = permute(temp, perm);
else
    gmdmproduct = diagsum(temp, d1, ndim1+d2);
end
 
% if ndims(tensor1) == 2 && d1 == 2 && ndims(tensor2) == 2 && d2 == 1
%     gmdmproduct = tensor1 * tensor2;
% else
%    
%     
%     if d1 == ndims(tensor1)+d2 && size(tensor2,d2) == 1
%         
%         perm = 1:ndims(temp);perm(d1) = []; perm = [perm, d1];
%         gmdmproduct = permute(temp,perm);
%     else
%         
%         gmdmproduct = diagsum(temp, d1, ndims(tensor1)+d2);
%     end
% end