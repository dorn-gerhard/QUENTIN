function [Current, Detail_Current, FL_I, FL_detail, Transmission, N] = CPT_current(G, bat, coupling, GF_grid, Mu, Beta, Bas, En)
%NOTE: is TT, TT_detail the transmission?

[B, coupl, G_hyb, Gam, TT, TT_detail, Trans_coh, Trans_incoh] = deal(cell(bat.Numel_Baths,1));

for bath_id = 1:bat.Numel_Baths
    
    
   
    s_func = (1-2*bat.Fermi(GF_grid.Points, Mu(bath_id,1), Beta(bath_id, 1), + 1));
    
    %NOTE: Initialzie Coupling Green functions and Bath Green functions
    
    if Bas.N_spin == 1
        B{bath_id,1} = Green(GF_grid, bat.f_GF_shifted(bath_id, Mu(bath_id,1)), [], s_func);
        coupl{bath_id} = coupling{bath_id,1} ;
        
    else
        
        if bat.Spin_conserved(bath_id) && bat.Spin_symmetric(bath_id)
        
        
            coupl{bath_id} = [coupling{bath_id,1}, zeros(Bas.N_site,bat.Numel_Orbitals(bath_id)); zeros(Bas.N_site, bat.Numel_Orbitals(bath_id)), coupling{bath_id,1}];
        
            B{bath_id,1} = Green(GF_grid, [bat.f_GF_shifted(bath_id,  Mu(bath_id,1)) , zeros(size(bat.GF{bath_id})); ...
                zeros(size(bat.GF{bath_id})), bat.f_GF_shifted(bath_id, Mu(bath_id,1)) ] , [], s_func);
        else
            coupl{bath_id} = coupling{bath_id,1};
            B{bath_id,1} = Green(GF_grid, bat.f_GF_shifted(bath_id, Mu(bath_id,1)));
        end
            
    end
    G_hyb{bath_id}.R = mmat(mmat(coupl{bath_id},B{bath_id}.R), coupl{bath_id}');
    G_hyb{bath_id}.K = mmat(mmat(coupl{bath_id},B{bath_id}.K), coupl{bath_id}');
    Gam{bath_id}.R = -1i/2 * (G_hyb{bath_id}.R - Green.f_dagger(G_hyb{bath_id}.R));
    Gam{bath_id}.K= -1i/2 * (G_hyb{bath_id}.K - Green.f_dagger(G_hyb{bath_id}.K));
end


FL_T = cell(bat.Numel_Baths,1);
for bath_id = 1:bat.Numel_Baths
    FL_T{bath_id} = zeros(size(G.R));
    Trans_coh{bath_id} = zeros(GF_grid.N,1);
    Trans_incoh{bath_id} = zeros(GF_grid.N,1);
    if bath_id == 1
        G_hybridization = Green(GF_grid, G_hyb{bath_id}.R, G_hyb{bath_id}.K);
    else
        G_hybridization.R = G_hybridization.R + G_hyb{bath_id}.R;
        G_hybridization.K = G_hybridization.K + G_hyb{bath_id}.K;
    end
end

G_full = Green.f_dyson(G, G_hybridization);

% Sigma.smaller = 1/2 * mmat(Green.f_inverse(G_full.R), mmat((G_full.K - (G_full.R - G_full.A)), Green.f_inverse(G_full.A)));
% Delta_Sigma = mmat(Green.f_inverse(G_full.R), mmat(G_full.R - G_full.A, Green.f_inverse(G_full.A)));
% Sigma.K = mmat(Green.f_inverse(G_full.R), mmat(G_full.K, Green.f_inverse(G_full.A)));
% %Sigma.comp = 1/2*(Sigma.K - Delta_Sigma);
%Test = Sigma.comp - Sigma.smaller;

%G lesser corresponds to <cdag c>
N = trace(Green.f_integrate_full(G_full.K - (G_full.R - G_full.A), GF_grid.Points)./(2i));
%N/(2*pi)
% plot_GF(diagsum(G_full.K - (G_full.R - G_full.A), 1,2), GF_grid.Points)

I_detail = cell(bat.Numel_Baths, 1);
FL_detail = cell(bat.Numel_Baths, 1);

if Bas.N_spin == 2 && ismember('S', En.Structure) %NOTE: spin conserved (in the bath and in the system and in the coupling
    %TODO: check coupling and spin conservation in bath
    numel_spin = 2; 
else
    numel_spin = 1;
end

I = zeros(bat.Numel_Baths, numel_spin);
FL_I = zeros(bat.Numel_Baths, numel_spin);

for bath_id = 1:bat.Numel_Baths
    
    %Variante 1 (\Delta * G)
    T = mmat(G_hyb{bath_id}.R, G_full.K) + mmat(G_hyb{bath_id}.K, Green.f_dagger(G_full.R)); 
    I_detail{bath_id} = Green.f_integrate_full(T, GF_grid.Points);
    
    %Variante 2 (Fisher Lee - landauer type current formula)
    other_baths = 1:bat.Numel_Baths;
    other_baths(bath_id) = [];
  
    for oth_bat = other_baths
        f_bath =  permute(bat.Fermi(GF_grid.Points, Mu(bath_id,1), Beta(bath_id, 1), 1), [2,3,1]);
        f_oth_bat = permute( bat.Fermi(GF_grid.Points, Mu(oth_bat,1), Beta(oth_bat, 1), 1),[2,3,1]);
        FL_T{bath_id} = FL_T{bath_id} + 4 * bsxfun(@times, f_bath - ...
            f_oth_bat , mmat(G_full.R , mmat(Gam{oth_bat}.R, mmat(G_full.A, Gam{bath_id}.R)))); ...
            
        % formula (6) from [Phys Rev Let. 1992, Meir Wingreen, Landauer Formula for the Current through an
        % Interacting Electron Region] 
        % Alternative for current
        TT{bath_id} =  + 1i/2*mmat((bsxfun(@times, f_bath, Gam{bath_id}.R) - bsxfun(@times, f_oth_bat, Gam{oth_bat}.R)), G_full.R - G_full.A) ... 
            + 1i/4 * mmat(Gam{bath_id}.R - Gam{oth_bat}.R, G_full.K - (G_full.R - G_full.A));
        
        Trans_coh{bath_id} = Trans_coh{bath_id} + diagsum(mmat(mmat(G_full.A, Gam{bath_id}.R), mmat(G_full.R, Gam{oth_bat}.R)),1,2);
        Trans_incoh{bath_id} = Trans_incoh{bath_id} + 4 * diagsum(bsxfun(@times, f_bath - ...
            f_oth_bat , mmat(G_full.R , mmat(Gam{oth_bat}.R, mmat(G_full.A, Gam{bath_id}.R)))),1,2) ./ ( f_bath(:) - f_oth_bat(:)/numel(other_baths));
    end
    FL_detail{bath_id} = Green.f_integrate_full(FL_T{bath_id}, GF_grid.Points);
    TT_detail{bath_id} = Green.f_integrate_full(TT{bath_id}, GF_grid.Points);
    
    %Transmission
   
    
    %separate by spin
    
    
    if ismember('S', En.Structure)
        for spin_id = 1:Bas.N_spin
            index = (1:Bas.N_site) + Bas.N_site * (spin_id - 1);
            I(bath_id, spin_id) = real(trace(I_detail{bath_id}(index, index)));
            FL_I(bath_id, spin_id) = real(trace(FL_detail{bath_id}(index, index)));
        end
    else
        I(bath_id, 1) = real(trace(I_detail{bath_id}));
        FL_I(bath_id, 1) = real(trace(FL_detail{bath_id}));
    end
end
Transmission = Trans_coh{1};

Current = I;
Detail_Current = I_detail;


% fisher lee type

%NOTE regarding transmission - the incoherent part cannot be reconstructed by division through the fermi functions
%since they tend to be zero in some ranges. -> rethink incoherent part, derive an explicit expression!
% plot_GF(Trans_coh{1}, GF_grid)
%plot_GF(Trans_incoh{1}, GF_grid)

%sum(abs(Trans_coh{2} - Trans_coh{1}))




 
 