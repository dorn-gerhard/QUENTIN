

%test for simple chain

t = 1
eps_0 = 1;
%expect semicircular DOS

V = t;
s = 0.;
tL = eps_0;
S = 1.*eye(length(eps_0));
omega_range = linspace(-6,6,1000);
NullPlus = 5;
epsilon = 10^-5;
nmax = 4000;


for omega_ind = 1:numel(omega_range)
    omega = omega_range(omega_ind);
    G_semi(omega_ind)  = Bath_New.sf_semi_chain_sancho( V,s,tL,S,omega,NullPlus,epsilon,nmax );
end

plot_GF(G_semi, omega_range)




%test for coupled chain
lead_coupling = 0.9;

t_leads = [1.3,            lead_coupling   ;...
           conj(lead_coupling), 0.3375]; %orbitals
eps_0 = [-2.1,   0   ; ... % shift in orbitals
          0  ,  -1.1];
%expect semicircular DOS

V = tril(t_leads);
s = zeros(size(t_leads));
tL = eps_0 + (t_leads-diag(diag(t_leads)));
S = 1.*eye(size(t_leads));
omega_range = linspace(-6,6,1000);
NullPlus = 5;
epsilon = 10^-5;
nmax = 4000;

clear G_semi
for omega_ind = 1:numel(omega_range)
    omega = omega_range(omega_ind);
    G_semi(:,:,omega_ind)  = Bath_New.sf_semi_chain_sancho( V,s,tL,S,omega,NullPlus,epsilon,nmax );
end

plot_GF(G_semi, omega_range)