 function result = ones_kron(x)
   
%    sx = size(x);
%    nx = prod(sx);
%    x = reshape(x, 1, nx);
%    
%    result = zeros(n^2, nx);
%    result(1:n+1:n^2,:) = x(ones(1,n),:);
%    result = reshape(result, [ n n sx ]);
%    result = permute(result, [1 3 2 4:2+length(sx) ]);
%    result = reshape(result, [ n*sx(1) n*sx(2) sx(3:end) ]);



% alternative:
%x = rand(2,3)
%n = 3;
[nx,ny] = size(x);
%row / columnindex of input data

B = repmat(x,nx,ny);

index1 = reshape(1:nx^2,nx,nx).';
index2 = reshape(1:ny^2,ny,ny).';

result = conj(B(index1, index2)) + B;


