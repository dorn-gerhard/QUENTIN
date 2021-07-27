
N_part = 4;
Spin = 0;

profile on
tic
for k = 1:10000
    Energ = En.g_En(N_part, Spin);
end
toc
profile viewer
% 7.139 sec

profile on
tic
for k = 1:10000
    Energ2 = En.Energies('NS', N_part, Spin);
end
toc
profile viewer
% 9.21 sec


