function S = quantum_relative_entropy(rho,sigma)

if ~all(size(rho) == size(sigma))
    error('dimensions do not match - S = +inf?')
end

%Tr(rho log2 rho) - Tr(rho log2 sigma)

[ev,ew] = eig(sigma);
tolerance = 10^-12;

if any(abs(imag(diag(ew)))>tolerance) || any(diag(ew)< -tolerance)
    error('density matrix is not valid')
end
ew = denoise(ew,tolerance);

l2ew = zeros(size(ew));
l2ew(ew> 10^-10) = log2(ew(ew> 10^-10));



%tr(rho log2 sigma)
entropy_rel = -trace(rho*ev*l2ew*ev');
S = entropy_rel  - entropy(rho);
S = denoise(S, tolerance);