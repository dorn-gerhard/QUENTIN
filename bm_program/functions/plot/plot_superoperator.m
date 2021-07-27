function plot_superoperator(K, Tab_cut, Energy_cut, Energy_cut2)

addpath /afs/itp.tugraz.at/user/dorn/Matlab/CPT/plot

if false
    addpath '/afs/itp.tugraz.at/user/dorn/Matlab/CPT/classes/'
    addpath '/afs/itp.tugraz.at/user/dorn/Matlab/CPT/classes/mmat'
    addpath '/afs/itp.tugraz.at/user/dorn/Matlab/CPT/classes/class_functions'
    addpath '/afs/itp.tugraz.at/user/dorn/Matlab/CPT/classes/Structure'
    addpath '/afs/itp.tugraz.at/user/dorn/Matlab/CPT/classes/multinv'
    addpath '/temp/dorn/run/benzene/'
    
    load('lindblad.mat', 'lindblad')
    struct2vars(lindblad)
end

if size(K,1) < 5000
   
    [ev, ew] = eig(full(K));
    ew = diag(ew);
   
else
    [ev, ew] = eigs(K, 1000, 10^-10);
    ew = diag(ew);
end


%find out range

x_min = [min(real(ew)), max(real(ew))];
y_min = [min(imag(ew)), max(imag(ew))];


%setup radius, when eigenvalues are counted twice

radius = min(diff(x_min), diff(y_min))/500;

%axis equal?

%find out how to count eigenvalues with a certain distance
[ew1,ew2] = meshgrid(ew);

distance = abs(ew1-ew2);
L = distance < radius;
frequency = sum(L,2) ;
frequency(frequency > 12) = 13;
spy(L)


markers = {'.','+', 'v', 's', 'p', 'h', '*', 'd', 'x', 'o','<', '>', '^'}; 
colors = lines(numel(markers));
plot_settings

hold off
for k = 1:numel(ew)
    plot(real(ew(k)), imag(ew(k)), 'marker', markers{frequency(k)}, 'linestyle', 'none', 'color', colors(frequency(k),:))
    hold on
end

for k = 1:numel(markers)
pp(k) = plot(-1,nan, markers{k}, 'markersize', 6, 'color', colors(k,:));
hold on
end
legend(pp, '1','2','3','4','5','6','7','8','9','10','11', '12', '> 12')
xlabel('real part')
ylabel('imaginary part')
title(['Spectrum of ', num2str(size(K,1)), ' x ', num2str(size(K,2)), ...
    ' operator \newlineradius for counting degenerate eigenvalues: ', num2str(radius)])
grid on

Tab_cut.Numel_Vector_Elements = sum(Tab_cut.Block_Size.^2)
index = Tab_cut.f_choi_index_permutation()


% examine time evolution
        
       

diag_index_shortened = Tab_cut.Diagonal_Index();
% a) random initial state:
initial_sigma = Sigma(Tab_cut.Energy, [], Tab_cut, Energy_cut, Energy_cut_2);
vec = initial_sigma.f_vector();

%b) maximum entangled state (aka Linksvakuum)
%max_ent_state = Tab_cut.f_vector(eye(Tab_cut.Numel_Datasets)/Tab_cut.Numel_Datasets);
%        initial_sigma = Sigma(En, max_ent_state , Tab_cut, Energy_cut, Energy_cut_2);
%        vec = initial_sigma.f_vector();


delta_t = 1;
n = 1000;

for k = 1:n
    disp(k)
    vec = expm(K*delta_t) * vec;
    trac(k) = sum(vec(diag_index_shortened));
end
subplot(2,1,1)
plot(1:1000, real(trac))
subplot(2,1,2) 
plot(1:1000, imag(trac))



% get representation

blocks = Tab_cut.f_vector_matrix2block_index(1:numel(vec));

N = Tab_cut.N_part(blocks);
S = Tab_cut.Spin(blocks);
if ismember('A', Tab_cut.Order)
    E = Tab_cut.Average_Energies(blocks);
elseif ismember('E', Tab_cut.Order)
    E = Tab_cut.Energies(blocks);
end

tag = arrayfun(@(x,y,z) ['N: ' num2str(x), ', S: ', num2str(y), ', E: ', num2str(z)], N,S,E, 'uniformoutput', false);
eigenvalue_tag = arrayfun(@(x) num2str(x), ew_sort, 'uniformoutput', false)
% examine eigenvectors

[en_sort, en_ix] = sort(E);

[ew_sort, ix] = sort(abs(ew));
low_eigenvectors = ev(en_ix,ix);
spy(abs(low_eigenvectors ) > 10^-6)
[energ, NN] = meshgrid(1:numel(abs(ew_sort)), 1:numel(N));
pp = pcolor(energ, NN, abs(low_eigenvectors));
set(pp, 'edgecolor', 'none');
set(gca,'ytick', 1:numel(ix))
set(gca,'yticklabel', tag)
set(gca, 'xtick', 1:numel(ix))
set(gca, 'xticklabel', eigenvalue_tag)

%lindblad = vars2struct('K', 'Tab_cut', 'En', 'Energy_cut', 'Energy_cut_2')
%save('lindblad.mat', 'lindblad')

index = Tab_cut.f_choi_index_permutation()



