 function y = eye_kron(n, x)
   
tmp = repmat({x},n,1);
y = blkdiag(tmp{:});