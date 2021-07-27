% inverse of a kronecker product, given A = B \otimes C
function C = kron_inv(A,B)
n = sqrt(length(A));

%n = 3;
%A = rand(n^2)



x = kron(ones(n), reshape(1:n^2,n,n));

y = kron(reshape(1:n^2,n,n), ones(n));
index = sub2ind([n^2,n^2], x,y);





%B = rand(3);

%C = rand(3);

%A = kron(B,C)


C = reshape(kron(B(:),eye(n^2))\A(index(:)), n,n);

