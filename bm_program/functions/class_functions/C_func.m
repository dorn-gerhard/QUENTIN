function C_val = C_func(om, tau, beta, mu_alpha, Gamma, site, spin)
% TODO: redesign input arguments

%mu_alpha and Gamma are defined for all leads (Gamma may be spin dependent)

%beta, mu_alpha, Gamma are all bath attributes - function should be part of Bath class
if nargin < 5 || isempty(Gamma), Gamma = 1; end

if nargin < 7 || isempty(spin), index = 1; spin = 0; end
if nargin < 6 || isempty(site), index = 1; site = 0; end

s_index = 0;%Bas.bin_index(site,spin);
if s_index == 0 || numel(Gamma) == 1; Gamma = Gamma(1); else Gamma = Gamma(index); end

if nargin < 4 || isempty(mu_alpha), mu_alpha = 1; end
if nargin < 3 || isempty(beta), beta = 30; end
if nargin < 2 || isempty(tau), tau = 1; end



% generalize to leads and use mu_vec and gamma_vec
% capable of using mu_alpha and Gamma vectors with two elements
C_val = Gamma(1).*(1/2* 1./(exp(tau*beta*(-tau*om-mu_alpha(1)))+1) + 1i/(2*pi) * real(psiz(1/2 + 1i/(2*pi) * beta*(om + tau*mu_alpha(1))))) + ...
       (numel(mu_alpha) > 1) * Gamma(end).*(1/2* 1./(exp(tau*beta*(-tau*om-mu_alpha(end)))+1) + 1i/(2*pi) * real(psiz(1/2 + 1i/(2*pi) * beta*(om + tau*mu_alpha(end)))));

