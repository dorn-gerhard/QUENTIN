% test to illustrate influence of dissipator function!

mu_alpha_range = linspace(-1.5,1.5,200);
Deff = linspace(-2, 2, 200);

dag = -1
c_spin = -1

for mu_index = 1:numel(mu_alpha_range)
    mu_index
    for Deff_index = 1:numel(Deff)
        for alpha = 1:2
            if dag == 1
                coupling_mu = conj(Coupling{alpha,Bat.Spin_fun(c_spin)});
                coupling_kappa = Coupling{alpha,Bat.Spin_fun(c_spin)};
            else
                coupling_mu = Coupling{alpha,Bat.Spin_fun(c_spin)};
                coupling_kappa = conj(Coupling{alpha,Bat.Spin_fun(c_spin)});
            end
            
            
            Bath3D_G1{alpha} = Bat.f_gamma(Deff(Deff_index), mu_alpha_range(mu_index), 30, alpha, c_spin, dag);
            Temp = coupling_mu*Bath3D_G1{alpha}* coupling_kappa';
            min_ew(mu_index, Deff_index, alpha) = min(eig(Temp)); %indexing: Er, El, mu, kappa, alpha
            max_ew(mu_index, Deff_index, alpha) = max(eig(Temp));
            
        end
    end
end


[xx,yy] = meshgrid(mu_alpha_range, Deff);

figure
plot_settings
subplot(2,1,1)

a = surf(xx, yy, sum(min_ew,3));
set(a, 'edgecolor', 'none')
view(2)
subplot(2,1,2)
b = surf(xx,yy, real(sum(max_ew,3)));
set(b, 'edgecolor', 'none')
view(2)



