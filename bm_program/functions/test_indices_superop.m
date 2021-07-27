

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
             RRR = reshape( bat.f_gamma(Emeshb(:) - Emesha(:), bat.Mu(alpha), bat.Beta(alpha), alpha, [],1),En.Numel_En,En.Numel_En);
             
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
    Am{alpha} = gmdmp(Q_tensor, 3,3,permute(V{alpha},[3,4,1,2]),3,4); %size(En, En, N_orbital)
end


%NOTE: this gives memory overflow, split into blocks!

GGG = zeros(En.Numel_En^2);
for alpha = 1:bat.Numel_Baths
    
    VI = sparse(Am{alpha}(:));
    WI = sparse(G{alpha}(:));
    
    
    Gamma = full(VI .* LLL.*(WI + WI.')/2 .* VI');
    Gamma2 = full(permute(conj(Am{alpha}), [3,4,1,2]).*(permute(G{alpha}, [3,4,1,2]) + G{alpha})/2 .* Am{alpha});
    % both are identical, LLL filter is missing!!!
    GGG = GGG + Gamma;
end


mat_ind = Tab.Matrix_Index;
index = sub2ind([En.Numel_En, En.Numel_En], mat_ind(:,1), mat_ind(:,2));
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




K2_ = K_full(index,index);

not_index = (1:En.Numel_En^2)';
not_index(index) = [];
K_reordered = K_full([index;not_index],[index;not_index]);

K_full_diag = K_full;
K_full_diag(not_index, not_index) = 0;
K_full_off = K_full - K_full_diag;


figure
subplot(1,2,1)
spy(K_full_diag,'b')
hold on
spy(K_full_off);

x=get(gca,'children');
set(x(2),'color','b')
set(x(1),'color',[0.5,0.5,0.5])

subplot(1,2,2)
spy(abs(K_full_diag([index;not_index],[index;not_index]))> 10^-5, 'b')
hold on
spy(abs(K_full_off([index;not_index],[index;not_index]))> 10^-5, 'k')
x=get(gca,'children');
set(x(2),'color','b')
set(x(1),'color',[0.5,0.5,0.5])

figure
subplot(2,2,1)
spy(abs(K))
subplot(2,2,2)
spy(abs(K)>10^-10)
subplot(2,2,3)
spy(abs(K)>10^-5)
subplot(2,2,4)
spy(abs(K)>10^-3)

%% Sigma
sigma_tolerance = 10^-11;
mes_flag = true,
[Sig, status] = get_sigma(K, En, Tab_cut, status, mes_flag, Energy_cut, Energy_cut_2, sigma_tolerance)
Sig.s_plot


% full Sigma
En = Hub.f_solve(energy.calc, energy.second_cut);

Tab2 = Table(En,'');
[Sig_full, status] = get_sigma(K_full, En, Tab2, status, mes_flag, Energy_cut, Energy_cut_2, sigma_tolerance);

%Frobeniusnorm:
sqrt(sum(sum(abs(Sig.f_matrix - Sig_full.f_matrix).^2)))

Sig_full.s_plot
[ev,ew] = eig(full(K_full));

L = abs(diag(ew)) <10^-10;
sum(L)
sig_vec = ev(:,L);
sigm = reshape(sig_vec,8,8);

%% energy selection
%restriction to those combinations which are relevant
% includes restriction for block diagonality of density matrix with respect to N_part and the creation annihilation
% combinations

% same as (Q_tensor == 1) .* ((sig-dag).*sig > 0) (creation or annihilation) but works also for unordered systems
LQ = cell(2,1);
for dag = 1:2
    
    LQ{dag,1} = false(size(Q_mat{1,1,1}));
    for N_site = 1:Bas.N_site
        LQ{dag,1} = logical(LQ{dag,1} + logical(abs(Q_mat{N_site,dag,1})));
    end
    
end




%restriction to smaller energies for quantum regression theorem!
s = 2; sig = 1; t = 1; tau = 1;

I_1 = En.Energies < energy.energy_cut
I_2 = En.Energies < energy.energy_cut_2
I1 = I_1 | I_2;
I = I_1;
I2 = true(1,numel(En.Energies));
LA = diag(I1) * LQ{s,sig} * diag(I1) * LQ{t,tau} * diag(I) * LQ{3-s,sig} * diag(I);

I_a = any(LA,2);
I_c = any(LA,1);

I_b = any(diag(I_a) * LQ{s,sig} * diag(I1),1);
I_d = any(diag(I) * LQ{3-s,sig} * diag(I_c),2);

[I_a, I_b', I_c', I_d]

L_i = diag(I_a) * LQ{s,sig} * diag(I_b);
L_j = diag(I_c) * LQ{s,sig} * diag(I_d);


%gmdmp(ones(1,1,En.Numel_En) .* gmdmp(LQ{s,sig},2,2,diag(I_a),2,1), 2,3,permute(LQ{t,tau},[3,1,4,2]),2,4)

% part A
LLA = LLA +kron(reshape(diag(I1)*LQ{s,sig}*diag(I1),[],1), ones(1,En.Numel_En^2)) .* ...
    kron(diag(I1)*LQ{t,tau}*diag(I), ones(En.Numel_En,En.Numel_En)) .* kron(ones(En.Numel_En^2,1), reshape(diag(I) * LQ{3-s, sig}'*diag(I), 1,[]));


[a,b,c,d] = ind2sub(numel(En.Energies) * ones(1,4), find(LLA));

% part B
LLB = kron(reshape((diag(I1)*LQ{3-s,sig}*diag(I2))',[],1), ones(1,En.Numel_En^2)) .* ...
    kron(reshape(diag(I2)*LQ{s,sig}*diag(I1),[],1), ones(1,En.Numel_En^2)) .* kron( diag(I1) * LQ{t, tau}*diag(I),ones(En.Numel_En,En.Numel_En));


[a,b,c,d] = ind2sub(numel(En.Energies) * ones(1,4), find(LLB));
LB = unique([a,b,d],'row');

% part C
LLC =  kron( diag(I1) * LQ{t, tau}*diag(I),ones(En.Numel_En,En.Numel_En)) .*kron(ones(En.Numel_En,1), kron((diag(I)*LQ{3-s,sig}*diag(I2))', ones(1,En.Numel_En)) ) ...
  .* kron(ones(En.Numel_En,1), kron((diag(I2)*LQ{s,sig}*diag(I)), ones(1,En.Numel_En)) ) ;
[a,b,c,d] = ind2sub(numel(En.Energies) * ones(1,4), find(LLC));
LC = unique([a,b,d],'row');



% calculate gamma

En.Energies(b) - En.Energies(a)

En.Energies(d) - En.Energies(c)




figure
spy(L_j)
hold on
spy(L_i)

%% constructed violation of positivity:

ew_delta = 0.0004;
[ev,EW] = eig(GGG);
%measure for non--Markovianity:
Non_Markov = -denoise(sum(EW(real(EW) < 0)), 10^-12)
%denoise(1/2 * sum(abs(diag(EW)) - diag(EW)), 10^-12)
L_ew = find(diag(real(EW)) == min(diag(real(EW))));

figure
plot_settings
symb = {'x', 'o', '<', 's', '*'};
for k = 1:5
    EW(L_ew, L_ew) = EW(L_ew, L_ew) - ew_delta;
    
    GGG_new = ev * EW * ev' ;
    
    
    
    subplot(2,1,1)
    ew = eig(GGG_new);
    plot(real(ew), imag(ew), 'x')
    title('spectrum of the Gamma matrix')
    axis equal
    grid on
    hold on
    
    
    
    
    Gamma_used = GGG_new;
    
    Temp = permute(sum(sum(reshape(Gamma_used.*kron(ones(En.Numel_En), eye(En.Numel_En)),8,8,8,8),3),1), [4,2,1,3]);
    GG = reshape(permute(reshape(Gamma_used, 8,8,8,8), [1,3,2,4]), 64,64) - 1/2*(kron(eye(En.Numel_En), Temp) + kron(transpose(Temp), eye(En.Numel_En))) ;
    ew = eig(GG);
    
    
    subplot(2,1,2)
    plot(real(ew), imag(ew), 'marker', symb{k},'linestyle', 'none')
    title('spectrum of the dissipator created from Gamma without the Hamiltonian part')
    axis equal
    grid on
    hold on
end


plot_pdf('spectrum_comparison_constructed_violation')


%% create adjacency matrix for spin up and spin down

E_cut = {-inf, energy.energy_cut};

E_cut = {-inf, -95};


b_blocks = Tab.Block_Index('A', E_cut);

b_npart = Tab.N_part(b_blocks);
b_spin = Tab.Spin(b_blocks);

GROUNDPATH = '/temp/dorn/run/bm_program/';


addpath([GROUNDPATH, 'functions/class_functions/KronProd/'])




b_blocks_cell = arrayfun(@(x) Tab.List_Index(x), b_blocks, 'uniformoutput', false);

% spin up
a_blocks_cell = arrayfun(@(x,y) Tab.List_Index('NSA', x+1,y+1, E_cut),b_npart, b_spin,'uniformoutput',false) ;

xx = cell2mat(cellfun(@(x,y) reshape(repmat(x,1, numel(y)), [],1), a_blocks_cell, b_blocks_cell, 'uniformoutput', false));

yy = cell2mat(cellfun(@(x,y) reshape(repmat(y.',numel(x),1),[],1), a_blocks_cell, b_blocks_cell, 'uniformoutput', false));

Ad_up = sparse(xx,yy,ones(size(xx)), En.Numel_En, En.Numel_En);
   spy(Ad_up+Ad_up') 
   
superindex_up = find(Ad_up+Ad_up');

   
   
% spin up
a_blocks_cell = arrayfun(@(x,y) Tab.List_Index('NSA', x+1,y-1, E_cut),b_npart, b_spin,'uniformoutput',false) ;

xx = cell2mat(cellfun(@(x,y) reshape(repmat(x,1, numel(y)), [],1), a_blocks_cell, b_blocks_cell, 'uniformoutput', false));

yy = cell2mat(cellfun(@(x,y) reshape(repmat(y.',numel(x),1),[],1), a_blocks_cell, b_blocks_cell, 'uniformoutput', false));

Ad_down = sparse(xx,yy,ones(size(xx)), En.Numel_En, En.Numel_En);

spy(Ad_down + Ad_down')

superindex_down = find(Ad_down+Ad_down');
numel(superindex_up)

superindex = [superindex_up; superindex_down];
numel(superindex)
%indices for one Gamma!

%Gamma space:
numel(superindex)^2
%compared to 4^6^4
% We do not want to calculate all those, as we do not need all those


% write down all four terms
%Energy terms:
Eall = En.Energies;
[Emeshx, Emeshy] = meshgrid(Eall, Eall);

[x,y] = meshgrid(1:En.Numel_En);
sig = sign(x-y); % lower triangle: s is negative, upper triangle: s is positive
%sig(1:10,1:10)
G = cell(bat.Numel_Baths,1);
for alpha = 1:bat.Numel_Baths
    for l = 1:bat.Numel_Orbitals(alpha)
        for k = 1:bat.Numel_Orbitals(alpha)
            G{alpha}(:,:,l,k) = bat.f_G_eval(sig.*(Emeshy - Emeshx));
            G{alpha}(:,:,l,k) = (1i*(G{alpha}(:,:,l,k) - conj(permute(G{alpha}(:,:,l,k),[1,2,4,3])))).*1./(exp(bat.Beta(alpha)*((Emeshy - Emeshx)- sig*bat.Mu(alpha))) + 1);
        end
    end
end

%check limits: should be between zero and 1/pi? at the moment -0.4 and 1/3????
%TODO: Check the bat.f_G_eval for large inputs that are larger than grid limits!!!
%GF may be wrong!!! not equal to 1

%Get the full Q matrices:

Q_mat = En.g_all_full_Qmat;
spy(Q_mat{1,2,2}) %site, dag (-1,1), spin (-1,1)
Q_tensor = zeros(numel(Eall), numel(Eall), basis.N_site); %(a,b,\mu)
for k = 1:basis.N_site
    for j= 1:2
        for l = 1:2
            Q_tensor(:,:,k) = full(Q_mat{k,1,1} + Q_mat{k,1,2} + Q_mat{k,2,1} + Q_mat{k,2,2});
        end
    end
end
%Test = cell2mat(Q_mat);

spy(Q_mat{1,1,1})
hold on 
spy(Q_mat{1,2,1},'c')
spy(Q_mat{1,1,2},  'r')
spy(Q_mat{1,2,2},  'g')

LQ = cell(2,2);
for dag = 1:2
    for spin = 1:2
        LQ{dag,spin} = false(size(Q_mat{1,1,1}));
        for N_site = 1:Bas.N_site
            LQ{dag,spin} = logical(LQ{dag,spin} + logical(Q_mat{N_site,dag,spin}));
        end
    end
end

s = 1; sig = 1; t = 2; tau = 1;
I1 = I_1 | I_2;
I = I_1;
I2 = true(1,64);
LA = diag(I1) * LQ{s,sig} * diag(I1) * LQ{t,tau} * diag(I) * LQ{3-s,sig} * diag(I);

I_a = any(LA,2);
I_c = any(LA,1);

I_b = any(diag(I_a) * LQ{s,sig} * diag(I1),1);
I_d = any(diag(I) * LQ{3-s,sig} * diag(I_c),2);

[I_a, I_b', I_c', I_d]

L_i = diag(I_a) * LQ{s,sig} * diag(I_b);
L_j = diag(I_c) * LQ{s,sig} * diag(I_d);


gmdmp(ones(1,1,64) .* gmdmp(LQ{s,sig},2,2,diag(I_a),2,1), 2,3,permute(LQ{t,tau},[3,1,4,2]),2,4)

% part A
LLA = kron(reshape(diag(I1)*LQ{s,sig}*diag(I1),[],1), ones(1,64^2)) .* ...
    kron(diag(I1)*LQ{t,tau}*diag(I), ones(64,64)) .* kron(ones(64^2,1), reshape(diag(I) * LQ{3-s, sig}'*diag(I), 1,[]));


[a,b,c,d] = ind2sub(numel(En.Energies) * ones(1,4), find(LLA));

% part B
LLB = kron(reshape((diag(I1)*LQ{3-s,sig}*diag(I2))',[],1), ones(1,64^2)) .* ...
    kron(reshape(diag(I2)*LQ{s,sig}*diag(I1),[],1), ones(1,64^2)) .* kron( diag(I1) * LQ{t, tau}*diag(I),ones(64,64));


[a,b,c,d] = ind2sub(numel(En.Energies) * ones(1,4), find(LLB));
LB = unique([a,b,d],'row');

% part C
LLC =  kron( diag(I1) * LQ{t, tau}*diag(I),ones(64,64)) .*kron(ones(64,1), kron((diag(I)*LQ{3-s,sig}*diag(I2))', ones(1,64)) ) ...
  .* kron(ones(64,1), kron((diag(I2)*LQ{s,sig}*diag(I)), ones(1,64)) ) ;
[a,b,c,d] = ind2sub(numel(En.Energies) * ones(1,4), find(LLC));
LC = unique([a,b,d],'row');



% calculate gamma

En.Energies(b) - En.Energies(a)

En.Energies(d) - En.Energies(c)




spy(L_i)
figure
spy(L_j)
hold on

% get coupling


%V{alpha}(site, orbital)
for alpha = 1:bat.Numel_Baths
    V{alpha}(:,:) = contact.coupling{alpha};
end

% permute, 
% diagsum(A, [2,3]) 2,3 =dimensions for trace
% gmdmp(A,i,N_i, B, j, N_j) A..tensor, B..tensor, i..dimension of tensor A, j..dimension of tensor B, N_i..number
% of dimensions of tensor A, N_j..number of dimensions of tensor B

for alpha = 1:bat.Numel_Baths
    Am{alpha} = gmdmp(Q_tensor, 3,3,permute(V{alpha},[3,4,1,2]),3,4); %size(En, En, N_orbital)
end

%

%%% NOW ALL OBJECTS for four dimensional Gamma tensor ready (for different alphas and orbitals!)

G, size(Am{1})


% perform operation below with blocks!!!
% find indexing scheme (i1, i2, w)



%NOTE: this gives memory overflow
for alpha = 1:bat.Numel_Baths
    
    VI = sparse(Am{alpha}(:));
    WI = sparse(G{alpha}(:));
    
    KR = KronProd({WI,WI.'});
    
    KS = kron_sparse(VI, VI.');
    
    Gamma = VI .* (WI + WI.')/2 .* VI';
    Gamma = permute(conj(Am{alpha}), [3,4,2,1]).*(permute(G{alpha}, [3,4,2,1]) + G{alpha})/2 .* Am{alpha};
end


%% Try indexing scheme
E_cut = {-inf,-3};
Tab.Average_Energies
b_blocks = Tab.Block_Index('A', E_cut);
I_1 = ismember(Tab.List_Index, Tab.List_Index(b_blocks));

E_cut = {-3,-1};
b_blocks = Tab.Block_Index('A', E_cut);
I_2 = ismember(Tab.List_Index, Tab.List_Index(b_blocks));

E_cut = {-1,inf};
b_blocks = Tab.Block_Index('A', E_cut);
I_3 = ismember(Tab.List_Index, Tab.List_Index(b_blocks));


L = logical(Q_tensor(:,:,1));


Tab_NS = Table(En, 'NS');


figure
plot_settings
ax = axes;

spy(L, 'ko', 4)
hold on
spy(L.*(I_1.*I_1'), 'b', 14)
spy(L.*(I_2.*I_2'), 'g', 14)
spy(L.*(I_3.*I_3'), 'r', 14)
spy(L.*(I_1.*I_2' + I_2.*I_1'), 'c', 14)
spy(L.*((I_1 + I_2).*I_3' + I_3 .*(I_1 + I_2)'), 'y',14)
pl = ax.Children;
clmap = colormap('lines');
%a = plot_colorbar();
%b = plot_colorbar2;
a(1,:) = [0    0.3922    1.0000];
a(2,:) = [0.2000    0.7961    0.7961];
a(3,:) = [ 1.0000    0.1987    0.0627];
b(1,:) = [ 0.8333    0.9500    1.0000];
b(2,:) = [1.0000    0.8987    0.7627]; 

%b(1,:) = (a(1,:) + a(25,:)) + [1,1,1] * 0.3;
%b(1,b(1,:) > 1) = 1;
col = {[0.7,0.7, 0.7], a(1,:), a(2,:), a(3,:), b(1,:), b(2,:)}

for k = 1:numel(pl)
    pl(k).Color = col{7-k};
end
Tab_NS.s_plot

LL = L;%zeros(size(Q_tensor(:,:,1)));
LL(~I_1, ~I_1) = false;
spy(LL)
E = true(64);
E(~I_1, ~I_1) = false
spy(E)
b_npart = Tab.N_part(b_blocks);
b_spin = Tab.Spin(b_blocks);


%% loop over restricted blocks of sigma

if ismember('A', Tab.Order)
    [energ_sort, energ_sort_index] = sort(Tab.Average_Energies);
else
    [energ_sort, energ_sort_index] = sort(Tab.Energies);
end
%[energ_sort, energ_sort_index]
energ_sort_index(energ_sort  > E_cut{2}) = []; %delete blocks that exceed Energy_cut

final_block_index = numel(energ_sort_index);
break_set = [20:10:numel(energ_sort_index), final_block_index]; %final_block_index; %
break_set = 10:1:final_block_index;

% use sorted average energy blocks which restrict sigma!!!


for n_block_index = 1:numel(energ_sort_index)
    n_block = energ_sort_index(n_block_index);
    n_part = Tab.N_part(n_block);
    spin = Tab.Spin(n_block);
    
    indices_a = Tab.List_Index(n_block);
    
    % other side:
    indices_c = [];
    for dag = -1:2:1
        for c_spin = -1:2:1
            indices_c = [indices_c; Tab.List_Index('NSA', n_part+dag, spin+dag*c_spin, E_cut)]; %check if sorted, should not be sorted
            
            
            %for s_dag = -1:2:1
            %    for s_spin = -1:2:1
            
            s_dag = 0;
            s_spin = 0;
            
            indices_b = Tab.List_Index('NSA', n_part+dag+s_dag, spin+dag*c_spin+s_dag * s_spin, E_cut)
                    
    
    
    
            indices_d = [];
    
            
        end 
    end
    
    
   
    
    
    
    
    
    
end
