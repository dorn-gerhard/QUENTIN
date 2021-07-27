function B = partial_trace(A, system_num, dim1)

% e.g. A     .. 6x6
% system_num .. 1 (system over which the trace is performed)
% dim1       ..  2
% result: 3x3 matrix

n = length(A);
dim2 = n/dim1;
if system_num == 1
    L = logical(kron( eye(dim1), true(dim2)));
   
    B = sum(reshape(A(L), dim2, dim2, dim1),3);
    
elseif system_num == 2
    L = logical(kron( true(dim1), eye(dim2)));
     B = squeeze(sum(reshape(A(L), dim1, dim2, dim1),2));
    
    
end



A = rand(6)
A = A+A'
