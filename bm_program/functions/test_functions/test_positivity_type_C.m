%% test positivity of Gamma for symmetrized master equation approach

% vary the following parameters:

%A) symmetry of applied potentials (in exec.bm line 235)
%B) symmetry of baths by replacing bath.GF (in exec.bm line 96)
%C) excitation behaviour 1:N, vs 2:2 (in setup)
%D) coupling (equal vs non-equal) replace contact.coupling (in exec.bm line 40)

%which system to take? 

%try spinless 3


ID = 'spinless_3_200121_0';
%A) do in exec.bm 
%B) one bath shifted one not, baths assymetric, 
%C) 2:3
%D) do in exec.bm


%A) do in exec.bm 
%B) one bath shifted one not, baths assymetric, 
%C) 1:3
%D) do in exec.bm


ID = 'spinless_3_200218_0';
% two orbitals in first bath!



setup_file = ['./data/setup_', ID, '.mat'];

k = 16;
    exec_bm(setup_file, num2str(k-1));



    S0 = Sigma(En,[],Tab);
    
    n = 400;
    time = linspace(0,10^4,n);
    for j = 1:n
    
    S1 = full(Tab.f_matrix(expm(K*time(j))*S0.f_vector()));
    T(:,j) = diag(S1);
    end
    
    

