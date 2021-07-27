function K_full = tensor_lindblad()

%% create list for all possible excitations for a given energy cut \mathcal{E}

%aim: create LLL
E_cut = {-inf, energy.energy_cut};

if ismember('A', Tab.Order)
    b_blocks = Tab.Block_Index('A', E_cut);
    e_tag = 'A';
else
    b_blocks = Tab.Block_Index('E', E_cut);
    e_tag = 'E';
end

b_npart = Tab.N_part(b_blocks);



% spinless
b_blocks_cell = arrayfun(@(x) Tab.List_Index(x), b_blocks, 'uniformoutput', false);

a_blocks_cell = arrayfun(@(x) Tab.List_Index(['N', e_tag], x+1, E_cut),b_npart, 'uniformoutput',false) ;

xx = cell2mat(cellfun(@(x,y) reshape(repmat(x,1, numel(y)), [],1), a_blocks_cell, b_blocks_cell, 'uniformoutput', false));

yy = cell2mat(cellfun(@(x,y) reshape(repmat(y.',numel(x),1),[],1), a_blocks_cell, b_blocks_cell, 'uniformoutput', false));

Ad_up = sparse(xx,yy,ones(size(xx)), En.Numel_En, En.Numel_En);
%   spy(Ad_up+Ad_up') 
   
superindex = find(Ad_up+Ad_up'); % finds all possible non zero excitations (creating or annihilating one particle
% corresponds to the superindex == find(Q_tensor(:,:,1))
T = false(En.Numel_En,En.Numel_En);
T(superindex) = true;

L_same_block = En.N_part - En.N_part' == 0;
L_A = L_same_block & abs(En.Energies - En.Energies')< 0.3; %not so good!

%LLL summarized the possible combinations of (ab,cd) eigenstates due to the following relation:
% K_{ab,cd}   <a| c^s |b > < b| \rho |c > <c| c^-s |d>
LLL = kron(T(:), ones(1,En.Numel_En^2)) .* kron(ones(En.Numel_En^2,1), reshape(T,1,[])) .* kron(ones(En.Numel_En), L_same_block) .*kron(L_same_block, ones(En.Numel_En));% kron(ones(1,En.Numel_En), kron(L_same_block, ones(En.Numel_En,1)));

% create all possible values for spectral function times fermi function the dagger (s) depends on whether one is in
% the lower or upper half (due to ordering according to particle number



% =================================== Q matrices ===============================================
Eall = En.Energies;

Q_mat = En.g_all_full_Qmat;
%spy(Q_mat{1,2,1}) %site, dag (-1,1), spin (-1,1)
Q_tensor = zeros(numel(Eall), numel(Eall), basis.N_site); %(a,b,\mu)
for k = 1:basis.N_site
    for j= 1:2
        for l = 1:2
            Q_tensor(:,:,k) = full(Q_mat{k,1,1} +  Q_mat{k,2,1} );
        end
    end
end


% get coupling

%V{alpha}(site, orbital)
for alpha = 1:bat.Numel_Baths
    V{alpha}(:,:) = contact.coupling{alpha};
end

for alpha = 1:bat.Numel_Baths
    Am{alpha} = gmdmp(Q_tensor, 3,3,V{alpha},1,2); %size(En, En, N_orbital)
end


% ==================================== DOS * fermi ==============================================







[Emeshb, Emesha] = meshgrid(Eall, Eall);

% < a | c_\kappa^-s | b > gamma(E_b - E_a)

[x,y] = meshgrid(1:En.Numel_En);
sig = sign(y-x); % lower triangle: s is negative, upper triangle: s is positive
%sig(1:10,1:10)
G = cell(bat.Numel_Baths,1);
S = cell(bat.Numel_Baths,1);

GGG = zeros(En.Numel_En^2);
HS = zeros(En.Numel_En);
for alpha = 1:bat.Numel_Baths
    S{alpha} = zeros(En.Numel_En, En.Numel_En, bat.Numel_Orbitals(alpha), bat.Numel_Orbitals(alpha));
    
    
    %G{alpha}(:,:,l,k) = bat.f_G_eval(sig.*(Emeshb - Emesha),alpha);
    %G{alpha}(:,:,l,k) = (1i*(G{alpha}(:,:,l,k) - conj(permute(G{alpha}(:,:,l,k),[1,2,4,3])))).*1./(exp(-bat.Beta(alpha)*((Emeshb - Emesha)- sig*bat.Mu(alpha))) + 1);
    % compare with Bat.f_gamma(obj, omega, mu_alpha, beta_alpha, alpha, c_spin, dag)
    %RRR = reshape( permute(bat.f_gamma(Emeshb(:) - Emesha(:), bat.Mu(alpha), bat.Beta(alpha), alpha, [],1),[3,1,2]),En.Numel_En,En.Numel_En,bat.Numel_Orbitals(alpha),[]);
    
    %G{alpha}(:,:,l,k)
    G{alpha} =  bat.f_gamma_mat((Emeshb - Emesha),alpha, -sig); %-sig since Eb-Ea
    %G{alpha}(:,:,l,k) - RRR
    
    % Do later:
    %S{alpha} =  bat.f_sigma_mat((Emeshb - Emesha),alpha, sig);
    % For now:
    %             SSS_plus = reshape(bat.f_sigma(Emeshb(:) - Emesha(:), bat.Mu(alpha), bat.Beta(alpha), alpha, [], 1),En.Numel_En,En.Numel_En);
    %
    %             S{alpha}(sig == 1) = SSS_plus(sig== 1);
    %             SSS_minus = reshape(bat.f_sigma(Emeshb(:) - Emesha(:), bat.Mu(alpha), bat.Beta(alpha), alpha, [], -1),En.Numel_En,En.Numel_En);
    %             S{alpha}(sig == -1) = SSS_minus(sig== -1);
    %
    
    
    VI_g = reshape(Am{alpha},En.Numel_En^2,1,[]);
    WI_g = reshape(G{alpha}, En.Numel_En^2, 1, bat.Numel_Orbitals(alpha),[]);
    
    
    Gamma = VI_g .* LLL.*(WI_g + permute(WI_g,[2,1,3,4]))/2 .* conj(permute(VI_g,[2,1,4,3]));
    GGG = GGG + sum(sum(Gamma,3),4);
    
    
    

    WI_s = zeros(En.Numel_En, En.Numel_En, bat.Numel_Orbitals(alpha), bat.Numel_Orbitals(alpha));
    VI = reshape(VI_g, En.Numel_En, En.Numel_En, bat.Numel_Orbitals(alpha)) ;
    %WI = -1i/2 * bat.f_sigma_mat(Ea-Eb, alpha,sig);
    L_plus = repmat(sig==1,1,1,bat.Numel_Orbitals(alpha), bat.Numel_Orbitals(alpha));
    L_minus = repmat(sig==-1,1,1,bat.Numel_Orbitals(alpha), bat.Numel_Orbitals(alpha));
    WI_s(L_plus) =  -1i/2 * permute(bat.f_sigma(Emesha(sig==1) - Emeshb(sig==1), bat.Mu(alpha), bat.Beta(alpha), alpha, [], 1),[3,4,1,2]);
    WI_s(L_minus) = -1i/2 * permute(bat.f_sigma(Emesha(sig==-1) - Emeshb(sig==-1), bat.Mu(alpha), bat.Beta(alpha), alpha, [], -1),[3,4,1,2]);
    for l = 1:bat.Numel_Orbitals(alpha)
        for p = 1:bat.Numel_Orbitals(alpha)
            temp = (VI(:,:,l).*WI_s(:,:,l,p)) * VI(:,:,p)';
            HS = HS + (temp + temp')/2; 
        end
    end
    

end



C_G = reshape(GGG, En.Numel_En, En.Numel_En, En.Numel_En, En.Numel_En);
CCC = reshape(permute(C_G, [1,3,2,4]), En.Numel_En^2, En.Numel_En^2);

K_full = -1i * (kron(eye(En.Numel_En), HS + diag(En.Energies)) - kron( transpose(HS) + diag(En.Energies), eye(En.Numel_En))) + ...
         CCC - 1/2*( kron(eye(En.Numel_En), transpose(diagsum(reshape(GGG,8,8,8,8),1,3))) + kron(diagsum(reshape(GGG,8,8,8,8),1,3), eye(En.Numel_En)));

     
     
mat_ind = Tab.Matrix_Index;
index = sub2ind([En.Numel_En, En.Numel_En], mat_ind(:,1), mat_ind(:,2));

K3 = K_full(index,index);

ew = eig(K_full);
positivity = - sum(ew(real(ew) > 0)) % should be zero
imaginary_part = sum(abs(imag(ew(real(ew)> 0)))) %should be zero

H_LS_2 = blkdiag(H_LS_calculated{1},H_LS_calculated{2},H_LS_calculated{3},H_LS_calculated{4})


% time evolution


r = Sigma(En, [], Tab);

tau = 0.00001;

C = reshape(permute(reshape(expm(K_full * tau), ones(1,4) * En.Numel_En), [1,3,2,4]), En.Numel_En^2, En.Numel_En^2);
ewc = eig(C);
ewg = eig(GGG);

compl_pos = (sum(abs(ewc))/En.Numel_En - 1) / tau;
markov = sum(abs(ewg(real(ewg) <0)));

markov - compl_pos * En.Numel_En/2


r2 = eig(full(Tab.f_matrix(expm(K2*1000) * r.f_vector)))



%% quantum regression theorem!!!

om = bat.Omega.Points;
i0plus = bat.Omega.Imag_Infinitesimal *1i;
    
G_qr = zeros(Bas.N_site * Bas.N_spin, Bas.N_site * Bas.N_spin, numel(om));

[V,E] = eig(full(K_full));
V_ = inv(V);

for ind_i = 1:Bas.N_site
    
    % {c_j^dag , \sigma } 1/(\omega + K + i0^+)  c_i
    
    % creator annihilator
    %c_i:
    c_i_all(:,ind_i) = Q_mat{ind_i,1}(:);
end
%c_j^\dag:




tic
for ind_k = 1:numel(om)
    %X2 = (eye(length(K_full))*(om(ind_k)+ i0plus) + K_full) \ c_i(:);
    %X2 = inv(eye(length(K_full))*(om(ind_k)+ i0plus) + K_full)
    X = V * diag(1./((om(ind_k)+i0plus) + diag(E))) * V_;
    
    index_set = 1:30;
    tic
    X3 = permute(gmdmp(V.*1./( (permute(om(index_set),[3,2,1]) + i0plus) + diag(E).'),2,3,V_,1,2),[1,3,2]);
    toc
    
    for ind_j = 1:Bas.N_site
        c_jdag = Q_mat{ind_j,2};
        anticomm = Sig.f_matrix * c_jdag + c_jdag*Sig.f_matrix;
        for ind_i = 1:Bas.N_site
            G_qr(ind_i,ind_j,ind_k) = full(transpose(anticomm(:)) * X* c_i_all(:,ind_i));
            
        end
        
    end
end
toc

