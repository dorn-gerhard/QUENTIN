function F = fidelity(rho,sigma)

tolerance = 10^-12;

sqrt_rho = rho^(1/2);

F = trace( (sqrt_rho*sigma*sqrt_rho)^(1/2))^2;
F = denoise(F, tolerance);