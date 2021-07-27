function C = choi(K)


[M,N] = size(K);
if (round(sqrt(M)) - sqrt(M)) > 10* eps ||  (M ~= N)
    error('Input Matrix needs to have dimension of n^2 x n^2')
end

n = sqrt(N);

C = reshape(permute(reshape(full(K),n,n,n,n), [1,3,2,4]), n^2, n^2);

