function index = choi_index(n)


x = kron(ones(n), reshape(1:n^2,n,n));

y = kron(reshape(1:n^2,n,n), ones(n));
index = sub2ind([n^2,n^2], x,y);
