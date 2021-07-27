function [K_full, K_small, I_full, Sig_full, status, index, index_off] = lindblad_tensor(En, Tab, bat, contact)



%TODO:
% 1. find an elegant way to calculate LLL (without energy cuts)
% 2. find an elegant way to deal with spin
% 3. think about energy cuts

%aim: create LLL - work around!!!
E_cut = {-inf, inf};

if ismember('A', Tab.Order)
    b_blocks = Tab.Block_Index('A', E_cut);
    e_tag = 'A';
elseif ismember('E', Tab.Order)
    b_blocks = Tab.Block_Index('E', E_cut);
    e_tag = 'E';
else
    b_blocks = Tab.Block_Index();
    e_tag = '';

end

b_npart = Tab.N_part(b_blocks);



% spinless
b_blocks_cell = arrayfun(@(x) Tab.List_Index(x), b_blocks, 'uniformoutput', false);
if ~isempty(e_tag)
    a_blocks_cell = arrayfun(@(x) Tab.List_Index(['N', e_tag], x+1, E_cut),b_npart, 'uniformoutput',false) ;
else
    a_blocks_cell = arrayfun(@(x) Tab.List_Index(['N'], x+1),b_npart, 'uniformoutput',false) ;
end

xx = cell2mat(cellfun(@(x,y) reshape(repmat(x,1, numel(y)), [],1), a_blocks_cell, b_blocks_cell, 'uniformoutput', false));

yy = cell2mat(cellfun(@(x,y) reshape(repmat(y.',numel(x),1),[],1), a_blocks_cell, b_blocks_cell, 'uniformoutput', false));

Ad_up = sparse(xx,yy,ones(size(xx)), En.Numel_En, En.Numel_En);
%   spy(Ad_up+Ad_up') 
   
superindex = find(Ad_up+Ad_up'); % finds all possible non zero excitations (creating or annihilating one particle
% corresponds to the superindex == find(Q_tensor(:,:,1))
T = false(En.Numel_En,En.Numel_En);
T(superindex) = true;

L_same_block = En.N_part - En.N_part' == 0;
%L_A = L_same_block & abs(En.Energies - En.Energies')< 0.3; %not so good!

%TODO: find a better strategy that also works for spin!!! 
%TODO: think about implementing an energy cut?
%LLL summarized the possible combinations of (ab,cd) eigenstates due to the following relation:
% K_{ab,cd}   <a| c^s |b > < b| \rho |c > <c| c^-s |d>
LLL = kron(T(:), ones(1,En.Numel_En^2)) .* kron(ones(En.Numel_En^2,1), reshape(T,1,[])) .* kron(ones(En.Numel_En), L_same_block) .*kron(L_same_block, ones(En.Numel_En));% kron(ones(1,En.Numel_En), kron(L_same_block, ones(En.Numel_En,1)));

% create all possible values for spectral function times fermi function the dagger (s) depends on whether one is in
% the lower or upper half (due to ordering according to particle number
% ==================================== DOS * fermi ==============================================
Eall = En.Energies;
[Emeshb, Emesha] = meshgrid(Eall, Eall);

% < a | c_\kappa^-s | b > gamma(E_b - E_a)

[x,y] = meshgrid(1:En.Numel_En);
sig = sign(x-y); % lower triangle: s is negative, upper triangle: s is positive
%sig(1:10,1:10)
G = cell(bat.Numel_Baths,1);
S = cell(bat.Numel_Baths,1);
for alpha = 1:bat.Numel_Baths
    S{alpha} = zeros(En.Numel_En, En.Numel_En, bat.Numel_Orbitals(alpha), bat.Numel_Orbitals(alpha));
    for l = 1:bat.Numel_Orbitals(alpha)
        for k = 1:bat.Numel_Orbitals(alpha)
            
            %G{alpha}(:,:,l,k) = bat.f_G_eval(sig.*(Emeshb - Emesha),alpha);
            %G{alpha}(:,:,l,k) = (1i*(G{alpha}(:,:,l,k) - conj(permute(G{alpha}(:,:,l,k),[1,2,4,3])))).*1./(exp(-bat.Beta(alpha)*((Emeshb - Emesha)- sig*bat.Mu(alpha))) + 1);
            % compare with Bat.f_gamma(obj, omega, mu_alpha, beta_alpha, alpha, c_spin, dag)
             %RRR = reshape( bat.f_gamma(Emeshb(:) - Emesha(:), bat.Mu(alpha), bat.Beta(alpha), alpha, [],1),En.Numel_En,En.Numel_En);
             
             %G{alpha}(:,:,l,k)
             G{alpha} =  bat.f_gamma_mat((Emeshb - Emesha),alpha, sig);
             %G{alpha}(:,:,l,k) - RRR
             
             % Do later:
             %S{alpha} =  bat.f_sigma_mat((Emeshb - Emesha),alpha, sig);
             % For now: 
             SSS_plus = reshape(bat.f_sigma(Emeshb(:) - Emesha(:), bat.Mu(alpha), bat.Beta(alpha), alpha, [], 1),En.Numel_En,En.Numel_En);
             
             S{alpha}(sig == 1) = SSS_plus(sig== 1);
             SSS_minus = reshape(bat.f_sigma(Emeshb(:) - Emesha(:), bat.Mu(alpha), bat.Beta(alpha), alpha, [], -1),En.Numel_En,En.Numel_En);
             S{alpha}(sig == -1) = SSS_minus(sig== -1);
             
        end
    end
end
% =================================== Q matrices ===============================================


Q_mat = En.g_all_full_Qmat;
%spy(Q_mat{1,2,1}) %site, dag (-1,1), spin (-1,1)
Q_tensor = zeros(numel(Eall), numel(Eall), En.Basis.N_site); %(a,b,\mu)
for k = 1:En.Basis.N_site
    for j= 1:2
        for l = 1:2
            Q_tensor(:,:,k) = full(Q_mat{k,1,1} +  Q_mat{k,2,1} );
        end
    end
end


% get coupling

%V{alpha}(site, orbital)
[V,Am] = deal(cell(1, bat.Numel_Baths));
for alpha = 1:bat.Numel_Baths
    V{alpha}(:,:) = contact.coupling{alpha};
end

for alpha = 1:bat.Numel_Baths
    Am{alpha} = gmdmp(Q_tensor, 3,3,permute(V{alpha},[3,4,1,2]),3,4); %size(En, En, N_orbital)
end


%NOTE: this gives memory overflow, split into blocks!

GGG = zeros(En.Numel_En^2);
Gamma = cell(bat.Numel_Baths,1); %NOTE: Split for spin current
for alpha = 1:bat.Numel_Baths
    
    VI = sparse(Am{alpha}(:));
    WI = sparse(G{alpha}(:));
    
    % Canonical Redfield Bloch!
    Gamma{alpha} = full(VI .* LLL.*(WI + WI.')/2 .* VI');
    %Gamma2 = full(permute(conj(Am{alpha}), [3,4,1,2]).*(permute(G{alpha}, [3,4,1,2]) + G{alpha})/2 .* Am{alpha});
    % both are identical, LLL filter is missing!!!
    GGG = GGG + Gamma{alpha};
end


mat_ind = Tab.Matrix_Index;
index = sub2ind([En.Numel_En, En.Numel_En], mat_ind(:,1), mat_ind(:,2));
index_off = (1:En.Numel_En^2)';
index_off(index) = [];

LL = false(En.Numel_En);
LL(index) = true;


[Eb, Ea] = meshgrid(En.Energies);
[y,x] = meshgrid(1:En.Numel_En);
sig = sign(x-y); % lower triangle: s is negative, upper triangle: s is positive


WI_ = cell(2,1);
HS_ = cell(2,1);
for alpha = 1:bat.Numel_Baths
    WI_{alpha} = zeros(En.Numel_En);
    VI = Am{alpha};
    %WI = -1i/2 * bat.f_sigma_mat(Ea-Eb, alpha,sig);
    WI_{alpha}(sig==1) =  -1i/2 * bat.f_sigma(Ea(sig==1) - Eb(sig==1), bat.Mu(alpha), bat.Beta(alpha), alpha, [], 1);
    WI_{alpha}(sig==-1) = -1i/2 * bat.f_sigma(Ea(sig==-1) - Eb(sig==-1), bat.Mu(alpha), bat.Beta(alpha), alpha, [], -1);
    
    HS_{alpha} = (VI.*WI_{alpha})*VI';
    HS_{alpha} = HS_{alpha} + HS_{alpha}';  
end

HS = LL.*(HS_{1} + HS_{2})/2;

C_G = reshape(GGG, En.Numel_En, En.Numel_En, En.Numel_En, En.Numel_En);
CCC = reshape(permute(C_G, [1,3,2,4]), En.Numel_En^2, En.Numel_En^2);

K_full = -1i * (kron(eye(En.Numel_En), HS + diag(En.Energies)) - kron( transpose(HS) + diag(En.Energies), eye(En.Numel_En))) + ...
     CCC - 1/2*( kron(eye(En.Numel_En), transpose(diagsum(reshape(GGG,En.Numel_En,En.Numel_En,En.Numel_En,En.Numel_En),1,3))) + kron(diagsum(reshape(GGG,En.Numel_En,En.Numel_En,En.Numel_En,En.Numel_En),1,3), eye(En.Numel_En)));

 
K_small = K_full(index,index);

Tab_full = Table(En,'N'); %TODO: Add spin
sigma_tolerance = 10^-10;
mes_flag = true;
status.eigenvalue_tolerance = 10^-11; 
[Sig_full, status] = get_sigma(K_small, En, Tab_full, status, mes_flag, inf, inf, sigma_tolerance);
%NOTE: Sigma full refers to the fact that no partial secular approximation
%or ERMEA was used, particle and spin conservation are assumed!

sig_tensor = repmat(sig,1,1,En.Numel_En, En.Numel_En);

I_full = zeros(bat.Numel_Baths,1); %TODO: Add spin
for alpha = 1:bat.Numel_Baths
    I_full(alpha)  = real(trace(transpose(diagsum(reshape(Gamma{alpha},En.Numel_En, En.Numel_En, En.Numel_En, En.Numel_En)  .* sig_tensor,1,3))*Sig_full.f_matrix));
end



