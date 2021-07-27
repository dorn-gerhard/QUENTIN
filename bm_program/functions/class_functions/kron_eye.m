 function result = kron_eye(x, n)
   
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
[yy,xx] = meshgrid(1:ny,1:nx);

additive = repmat(1:n,nx*ny,1);

x_new = repmat(xx(:)-1,1,n)*n + additive;
y_new = repmat(yy(:)-1,1,n)*n + additive;
data_new = repmat(x(:),1,n);

result = sparse(x_new,y_new,data_new, n*nx, n*ny);

result = full(result);


