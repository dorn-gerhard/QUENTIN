function [A,B,C,D,E,F] = func_index_transform(N_part, N_site, Spin, tau, c_site, c_spin, block)
% creates a unique index for the given input parameters and transforms back
% used to identify which blocks are effectively used - shall save time
calc_index = true; 
if nargin < 3, calc_index = false; end

%Test: 
% 
% N_site = 2;
% A = 2; B = 2; C = 1; D = 2; E = -1; F = 3;
% X = [A,B,C,D,E,F];

N_A = 2*N_site + 1;
N_B = 2*N_site + 1;
N_C = 2;
N_D = N_site;
N_E = 2;
N_F = 5; % #5 diag term


if calc_index == true
    
    
    A = N_part;
    B = Spin;
    C = tau;
    D = c_site;
    E = c_spin;
    F = block;
    
    A = 1 + (F-1)  +  N_F * (E+1)/2  +  N_F*N_E * (D-1)  +  N_F*N_E*N_D * (C+1)/2  +  N_F*N_E*N_D*N_C * (B + N_site)  +  N_F*N_E*N_D*N_C*N_B * A;
    
else
    
    index = N_part;
    A = floor((index-1)/(N_B*N_C*N_D*N_E*N_F));
    B = floor(rem((index-1), N_B*N_C*N_D*N_E*N_F) / (N_C*N_D*N_E*N_F))  - N_site;
    C = floor(rem((index-1), N_C*N_D*N_E*N_F) / (N_D*N_E*N_F)) *2 -1;
    D = floor(rem((index-1), N_D*N_E*N_F) / (N_E*N_F)) +1;
    E = floor(rem((index-1), N_E*N_F) / N_F) *2-1;
    F = rem((index-1), N_F) + 1;
end
% Test
% [A,B,C,D,E,F] - X
