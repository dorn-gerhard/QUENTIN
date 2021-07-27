function S = entropy(rho)

[ev,ew] = eig(rho);
tolerance = 10^-12;

if any(abs(imag(diag(ew)))>tolerance) || any(diag(ew)< -tolerance)
    error('density matrix is not valid')
end
ew = denoise(ew,tolerance);

%ew(ew> 10^-10) = log2(ew(ew> 10^-10));
%S = -trace(rho*ev*ew*ev');

%work via eigenvalues:


ew2 = ew(ew> tolerance);

S = -sum( ew2 .* log2(ew2));