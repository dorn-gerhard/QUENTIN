function plot_results(file_name_result, print_now)

if nargin < 2 || isempty(print_now) , print_now = false; end

load(file_name_result)


struct2vars(data_combined)
struct2vars(parameters)
%% create ui_figure
fig = figure;
plot_settings;
tabgrp = uitabgroup(fig);

% visualize
%Hub_spectrum(parameters.basis, parameters.system)
%load(energy.file{1})
%En.s_spectrum([],[],1)
%plot_parameters(parameters)

%% Current tab - display all currents (2D or 1D)
cur_tab = uitab(tabgrp, 'Title', 'Current', 'backgroundcolor', 'w');
axes('parent', cur_tab)
tabgrp.SelectedTab = cur_tab;

h = 6.62607015 * 10^-34;
e = 1.602176634 * 10^-19;
e_factor =  2*e^2/h*10^6; %µA

I       = II         * e_factor;
I_cpt   = II_cpt     * e_factor;
I_bmcpt = II_bmcpt   * e_factor;


if isfield(contact, 'couple_bias_gate') && contact.couple_bias_gate
    num_Vg = 1;
else
    num_Vg = numel(Vg_vec); % both cases (single Vg and multiple Vg s)
end
num_Vb = numel(Vb_vec);



VB = repmat(Vb_vec, num_Vg, 1);
VG = repmat(Vg_vec.', 1, num_Vb);



Vb_delta = VB(2)- VB(1);
VB_delta = ones(size(VB))*Vb_delta;

if size(VB,1) > 1 && size(VB,2) > 1
    %3D plot
    
    
    subplot(3,1,1)
    plot_settings
    
    surf(VB, VG, I(:,:,2),'edgecolor', 'none');
    view(2)
    
    plot_supertitle(['$\verb|', condor.setup_id, '|$'], cur_tab)
    title('BM current')
    xlabel('V_b'), ylabel('V_g')
    hold on
    grid on
    plot_colorbar
    caxis([-max(abs(I(:))), max(abs(I(:)))])
    
    subplot(3,1,2)
    
    surf(VB, VG,I_cpt(:,:,2),'edgecolor', 'none');
    view(2)
    
    title('CPT current')
    xlabel('V_b'), ylabel('V_g')
    hold on
    grid on
    plot_colorbar
    caxis([-max(abs(I_cpt(:))), max(abs(I_cpt(:)))])
    
    
    subplot(3,1,3)
    xweight = cellfun(@(x) sum(x.Tab.XWeight), Sigg);
    surf(VB, VG, xweight,'edgecolor', 'none');
    view(2)
    title('Off diagonal Weight')
    xlabel('V_b'), ylabel('V_g')
    hold on
    grid on
    plot_colorbar
    caxis([-max(xweight(:)), max(xweight(:))])
    
    
    
    
    
else
    
    
    subplot(3,1,1)
    plot_settings
    
    plot(VB,-I(:,:,1),'-x');
    hold on
    plot(VB,I(:,:,2)','-x');
    plot(VB,-I(:,:,3),'-x');
    lab_bm = plot(VB,I(:,:,4),'-x');
    
    plot_supertitle(['$\verb|', condor.setup_id, '|$'], cur_tab)
    title('BM current')
    xlabel('V_b'), ylabel('I | µA')
    hold on
    grid on
    legend(lab_bm, 'bm current', 'location', 'nw')
    
    subplot(3,1,2)
    
    plot(VB,-I_cpt(:,:,1),'-o');
    hold on
    plot(VB,I_cpt(:,:,2),'-o');
    plot(VB,-I_cpt(:,:,3),'-o');
    lab_cpt = plot(VB,I_cpt(:,:,4),'-o');
    
    
    title('CPT+BMCPT Current')
    xlabel('V_b'), ylabel('I | µA')
    grid on
    
    
    
    plot(VB,-I_bmcpt(:,:,1),'--*');
    hold on
    plot(VB,I_bmcpt(:,:,2),'--*');
    plot(VB,-I_bmcpt(:,:,3),'--*');
    lab_bmcpt = plot(VB,I_bmcpt(:,:,4),'--*');
    
    legend([lab_cpt, lab_bmcpt], 'cpt current', 'bmcpt current', 'location', 'nw')


    subplot(3,1,3)
    L = false(size(Sigg));
    L(~cellfun(@(x) isempty(x), Sigg)) = true;
    xweight = nan(size(Sigg));
    xweight(L) = cellfun(@(x) sum(x.Tab.XWeight), Sigg(L));
    plot(VB, xweight, '-.o');
    title('Off diagonal Weight')
    ylabel('Off diagonal Weight'), xlabel('V_b')
    grid on
end
if print_now
    plot_pdf(condor.setup_id)
end

%% plot parameters
param_tab = uitab(tabgrp, 'Title', 'Parameters', 'backgroundcolor', 'w');
tabgrp.SelectedTab = param_tab;
text_param = plot_parameters(parameters, '\n');
annotation(param_tab, 'textbox',[0.05 0.05 1 0.95], 'String', text_param,'FitBoxToText','on', 'fontsize', 8);

if print_now
    plot_pdf(condor.setup_id)
end


%% Conductivity
conduc_tab = uitab(tabgrp, 'Title', 'Conductivity', 'backgroundcolor', 'w');
axes('parent', conduc_tab)
tabgrp.SelectedTab = conduc_tab;


% use higher numeric differentiation method than simple finite differences
% NOTE: use http://web.media.mit.edu/~crtaylor/calculator.html to get the coefficients
dif_I_bmcpt = finite_derivative(I_bmcpt(:,:,2),2,Vb_delta);

%dif_I_bmcpt = diff(I_bmcpt(:,:,2), [],2)./Vb_delta;
dif_I_cpt = finite_derivative(I_cpt(:,:,2), 2, Vb_delta);
dif_I_bm = finite_derivative(I(:,:,2), 2, Vb_delta);


if size(VB,1) > 1 && size(VB,2) > 1
    subplot(3,1,1)
    
    
    surf(VB, VG, dif_I_bm,'edgecolor', 'none');
    view(2)
    
    plot_supertitle(['$\verb|', condor.setup_id, '|$'], conduc_tab)
    title('BM differential conductance')
    xlabel('V_b'), ylabel('V_g')
    hold on
    grid on
    plot_colorbar
    caxis([-max(abs(dif_I_bm(:))), max(abs(dif_I_bm(:)))])
    
    subplot(3,1,2)
    
    surf(VB, VG, dif_I_cpt,'edgecolor', 'none');
    view(2)
    
    title('CPT differential conductance')
    xlabel('V_b'), ylabel('V_g')
    hold on
    grid on
    plot_colorbar
    caxis([-max(abs(dif_I_cpt(:))), max(abs(dif_I_cpt(:)))])
    
    
    subplot(3,1,3)
    
    surf(VB, VG, dif_I_bmcpt,'edgecolor', 'none');
    view(2)
    title('BMCPT differential conductance')
    xlabel('V_b'), ylabel('V_g')
    hold on
    grid on
    plot_colorbar
    caxis([-max(dif_I_bmcpt(:)), max(dif_I_bmcpt(:))])
    
    
    
    
   
    
    

else
    subplot(3,1,1)
    
    
    lab_bm = plot(VB,dif_I_bm,'-x');
    
    
    
    plot_supertitle(['$\verb|', condor.setup_id, '|$'], conduc_tab)
    title('BM differential conductance')
    xlabel('V_b'), ylabel('G')
    hold on
    grid on
    legend(lab_bm, 'bm conductance', 'location', 'nw')
    

    subplot(3,1,2)
    
    lab_cpt = plot(VB,dif_I_cpt,'-o');
    
    title('CPT + BMCPT differential conductance')
    xlabel('V_b'), ylabel('G')
    grid on
    hold on
    lab_bmcpt = plot(VB,dif_I_bmcpt,'--*');
    
    legend([lab_cpt, lab_bmcpt], 'cpt conductance', 'bmcpt conductance', 'location', 'ne')
    
    
    subplot(3,1,3)
    
    pc = plot(VB,NN, '-x');
    plot_supertitle(['$\verb|', condor.setup_id, '|$'])
    title('Particle sector')
    %plot_colorbar_integer;
    %caxis([-1, 1]*30)
    xlabel('V_b'), ylabel('N')
    grid on
    

end

if print_now
     plot_pdf(condor.setup_id)
end



%% analysis of data quality

data_qual_tab = uitab(tabgrp, 'Title', 'Data Quality', 'backgroundcolor', 'w');
axes('parent', data_qual_tab)
tabgrp.SelectedTab = data_qual_tab;



if size(VB,1) > 1 && size(VB,2) > 1
    %3D plot
    subplot(3,2,[1,2])
    plot_supertitle(['$\verb|', condor.setup_id, '|$'], data_qual_tab)
    surf(VB, VG, log10(column_sum), 'edgecolor', 'none')
    view(2)
    
    title('column sum (logarithmic)')
    xlabel('V_b'), ylabel('V_g')
    hold on
    grid on
    plot_colorbar
    caxis([-max(abs(log10(column_sum(:)))), max(abs(log10(column_sum(:))))])
    
    subplot(3,2,[3,4])
    surf(VB, VG, log10(minimal_value), 'edgecolor', 'none')
    view(2)
    
    title('minimal eigenvalue (logarithmic)')
    xlabel('V_b'), ylabel('V_g')
    hold on
    grid on
    plot_colorbar
    caxis([-max(abs(log10(minimal_value(:)))), max(abs(log10(minimal_value(:))))])
    
    subplot(3,2,5)
    surf(VB, VG, multiple_solutions, 'edgecolor', 'none')
    view(2)
    
    title('multiple solutions')
    xlabel('V_b'), ylabel('V_g')
    hold on
    grid on
    plot_colorbar_integer
    caxis([-max(abs(multiple_solutions(:))), max(abs(multiple_solutions(:)))])
    
    subplot(3,2,6)
    surf(VB, VG, positive_eigenvalues, 'edgecolor', 'none')
    view(2)
    
    title('positive eigenvalues')
    xlabel('V_b'), ylabel('V_g')
    hold on
    grid on
    plot_colorbar
    caxis([-max(abs(positive_eigenvalues(:))), max(abs(positive_eigenvalues(:)))])
    


else
    subplot(3,1,1)
    plot_supertitle(['$\verb|', condor.setup_id, '|$'], data_qual_tab)
    
    plot(VB,log10(column_sum), '-x')
    title(' Columns sum (logarithmic)')
    xlabel('V_b'), ylabel('column\_sum (log10)')
    grid on
    
    
    subplot(3,1,2)
    
    leg_min_ev = plot(VB,log10(abs(minimal_value)), '-x');
    title('minimal eigenvalue and tolerance (log10)')
    hold on
    %plot_colorbar;
    xlabel('V_b'), ylabel('minimal eigenvalue')
    
    leg_lim_ev = plot(VB,log10(abs(eigenvalue_tolerance)), '-o');
    grid on
    legend([leg_min_ev,leg_lim_ev], 'minimal eigenvalue', 'eigenvalue tolerance', 'location', 'NW')
    
    subplot(3,1,3)
    
    yyaxis left
    lab_multi = plot(VB,multiple_solutions, '-x');
    title('Positive Eigenvalues | Multiple Solutions')
    %plot_colorbar_integer(8);
    hold on
    ylabel('Multiple solutions')
    ylim([0,max(multiple_solutions)+1])
    
    yyaxis right
    lab_pos_ev = plot(VB,positive_eigenvalues, '-o');
    
    xlabel('V_b'), ylabel('Positive eigenvalues')
    grid on
    legend([lab_multi, lab_pos_ev], 'multiple solutions', 'positive eigenvalues')
end


if print_now
    plot_pdf(condor.setup_id)
end

%% third evaluation


data_occu_tab = uitab(tabgrp, 'Title', 'Occupation', 'backgroundcolor', 'w');
axes('parent', data_occu_tab)
tabgrp.SelectedTab = data_occu_tab;

if size(VB,1) > 1 && size(VB,2) > 1
    %3D plot
    %TODO: think about it - 3D scatter plot? choose a line?
else
    
    
    plot_supertitle(['$\verb|', condor.setup_id, '|$'], data_occu_tab)
    subplot(1,9,[1,8])
    plot_settings
    plot(VB, data_combined.energy_cut, 'r')
    hold on
    %legend('maximal_calc_energy')
    grid on
    
    xlabel('VB')
    ylabel('occupied eigenstates')
    title(' occupation of eigenstates in steady state dependend on applied voltage')
    %x = colormap(parula(100));
    
    %colormap for 10^-1 range
    %x = colormap(gca, hsv(130))  ;
    x = flipud(colorcet('L9', 'N', 130));
    x = x(31:130,:);
    %colormap(gca, x)
    
    %grey colormap
    %start_color = [0.48,0.77,0.80];
    %y = flipud([linspace(start_color(1),1,100)', linspace(start_color(2),1,100)', linspace(start_color(3),1,100)']);
    %y = flipud(hot(100));
    y = flipud(colorcet('L8', 'N', 130));
    cmap = [y(1:100,:);x];
    colormap(gca,cmap)
    
    
    % gives all probalities a number between zero and 20
    min_val = 10^-6;
    linlogscale_fun = @(x) 0.00001 + (x> min_val) .*( (10+ (1+log10(x)) * 10/ (-log10(min_val)-1)) .* (log10(x) < -1) + (10*x+9) .* (log10(x) >= -1));
    
    for k = 1:numel(VB)
        %Average vs Energies
        if ~isempty(Sigg{k})
            T = data_combined.Sigg{k}.s_dom_blocks;
            if ismember('A', data_combined.Sigg{k}.Block_structure)
                EE = data_combined.Sigg{k}.s_dom_blocks.Average_Energies;
            else
                EE = data_combined.Sigg{k}.s_dom_blocks.Energies;
            end
            
            %load(data_combined.parameters.energy.file{k}, 'En')
            if isfield(contact, 'couple_bias_gate') && contact.couple_bias_gate
                if isfile(energy.file{Vg_vec == Vg_vec(k)})
                    load(energy.file{Vg_vec == Vg_vec(k)})
                else
                    En = get_En(setup_file,Vg_vec(k));
                end
                %load(['./data/energy_Vg_', num2str(Vg_vec(k)), '.mat'])
            else
                if k == 1
                    if isfile(energy.file{1})
                        load(energy.file{1})
                    else
                        En = get_En(setup_file, Vg_vec(1));
                    end
                end
                %load(['./data/energy_Vg_', num2str(Vg_vec(1)), '.mat'])
            end
            energies = En.Energies('E', {-inf, data_combined.energy_cut(k)});
            
            %plot(repmat(VB(k), numel(energies), 1), energies, 's', 'markersize', 5, 'color', [1,1,1]*0.70)

%             for j = 1:numel(energies)
%                 plot(VB(k), energies(j), 's', 'markersize', 5, 'color', [1,1,1]*0.70)
%             end
            % get degeneracies for EE
            round_E = round(EE*10^6)/10^6;
            [unique_E, ~, iy] = unique(round_E);
            a = hist(round_E, unique_E);            
            degeneracy = a(iy)';

            for j = 1:length(T.Index)
                
                WW = T.Weight(j);
                if T.Weight(j) > 1
                    warning ('probability > 1!')
                else
                    %previously used color: x(ceil(WW*length(x)),:)8+log10(WW)
                    plot(VB(k), EE(j), 'o', 'markersize', 4, 'markerfacecolor', cmap(ceil(linlogscale_fun(WW)*10),:), 'color', cmap(ceil(linlogscale_fun(WW)*10),:))
                    if degeneracy(j) > 1
                       % plot(VB(k), EE(j), '.k')
                    end
                end
                hold on
            end
            
        end
    end
    subplot(1,9,9)
    image(permute(cmap, [1,3,2]))
    
    factor = -log10(min_val);
    set(gca,'ydir', 'normal')
    set(gca,'Ytick', [linspace(0.5,100.5, factor), linspace(111,200,9)])
    set(gca, 'Yticklabel', cellfun(@(x) num2str(x),num2cell([10.^[-factor:-1], 0.2:0.1:1]), 'uniformoutput', false))
    set(gca, 'xtick', [])
    set(gca,'YAxisLocation', 'right')
    
    
end

if print_now 
    plot_pdf(condor.setup_id)
end
%% show setup + spectrum
system_tab = plot_setup(parameters.condor.setup_id, print_now, tabgrp);
axes('parent', system_tab)
subplot(2,1,2)

%% show spectrum

% data_spec_tab = uitab(tabgrp, 'Title', 'Spectrum', 'backgroundcolor', 'w');
% axes('parent', data_spec_tab)
% tabgrp.SelectedTab = data_spec_tab;

En.s_spectrum([],[],1, gca)
plot_supertitle(['$\verb|', condor.setup_id, '|$'], system_tab)

if print_now 
    plot_pdf(condor.setup_id)
end

%% Define Research Answer to Research Question

if ~isfield(data_combined.parameters.condor, 'research_question')
    data_combined.parameters.condor.research_question = '';
end
if ~isfield(data_combined.parameters.condor, 'research_answer')
    data_combined.parameters.condor.research_answer = '';
end
prompt = {'Research question:','Research answer:'};
dlgtitle = 'Documentation';
dims = [1 200];
definput = {data_combined.parameters.condor.research_question, data_combined.parameters.condor.research_answer};
options.Resize='off';
options.WindowStyle='normal';
answer = inputdlg(prompt,dlgtitle,dims,definput, options);
if ~isempty(answer)
    data_combined.parameters.condor.research_question = answer{1};
    data_combined.parameters.condor.research_answer = answer{2};
    save(file_name_result, 'data_combined')
end


end
%{
%%  generate table with all relevant eigenstates
last_successful_index = find(data_combined.no_solution == 0,1, 'last');
full_table = Sigg{last_successful_index}.s_dom_blocks;
old_ind = Sigg{last_successful_index}.Tab.Old_Indices(Sigg{last_successful_index}.f_diagonal_weight> 10^-5);
Sigma_qualification = zeros(numel(Sigg), numel(old_ind), 4); %weight, neg.eigenvalues, neg diag, factor of entanglement
for k = 1:last_successful_index
    for l = 1:numel(old_ind)
        Sigma_block = Sigg{k}.Values{Sigg{k}.Tab.Block_Index('O', old_ind(l))};
        ev = eig(full(Sigma_block));
        corr_meas = max(sum(abs(Sigma_block - diag(diag(Sigma_block))),1)); %is hermitian
        Sigma_qualification(k,l,:) = [sum(abs(diag(Sigma_block))), sum(real(ev) < 0), sum(diag(Sigma_block) < 0), corr_meas]; 
    end
end


figure
plot_settings
subplot(3,1,1)
plot(Vb_vec, Sigma_qualification(:,:,1)')
if ismember('E', Sigg{1}.Tab.Categ.f_shortnames)
    legend(num2str(Sigg{last_successful_index}.Tab.Energies('O', old_ind)))
    title('Occupation of density matrix versus applied voltage V_B for different Eigenstates (E)')
else
    legend(num2str(Sigg{last_successful_index}.Tab.Average_Energies('O', old_ind)))
    title('Occupation of density matrix versus applied voltage V_B for different Averaged Eigenstates (A)')
end
xlabel('V_B')
ylabel('Occupation')
plot_supertitle(['$\verb|', condor.setup_id, '|$'])

grid on

subplot(3,1,2)
plot(Vb_vec, Sigma_qualification(:,:,4)')
if ismember('E', Sigg{1}.Tab.Categ.f_shortnames)
    legend(num2str(Sigg{last_successful_index}.Tab.Energies('O', old_ind)))
    title('Correlation of density matrix versus applied voltage V_B for different Eigenstates (E)')
else
    legend(num2str(Sigg{last_successful_index}.Tab.Average_Energies('O', old_ind)))
    title('Correlation of density matrix versus applied voltage V_B for different Averaged Eigenstates (A)')
end
grid on
xlabel('V_B')
ylabel('Correlation')


subplot(3,1,3)
plot(Vb_vec, sum(Sigma_qualification(:,:,2)'))
hold on
plot(Vb_vec, sum(Sigma_qualification(:,:,3)'))
legend('negative eigenvalue in density matrix', 'negative entries on diagonal')
grid on
xlabel('V_B')
ylabel('Bad density matrix')


if print_now 
    plot_pdf(condor.setup_id)
end


%% check special point

k = 12
ID = data_combined.parameters.condor.setup_id
load(['./log_files/K_', ID, '_1x', num2str(k), '.mat'])
struct2vars(K_data)
figure
spy(K_data.K)
[ev,ew] = eig(full(K_data.K));
ew = diag(ew);
L = abs(ew) < 10^-10;
ev(:,L)
figure
plot(real(ew), imag(ew), 'x')
 [ev(:,L), K_data.Tab_cut.Average_Energies, K_data.Tab_cut.N_part, K_data.Tab_cut.Spin, K_data.Tab_cut.Block_Size]

 
 [Sig2, status] = get_sigma(K, Tab_cut, status, false, Energy_cut, Energy_cut_2, 10^-10)
 Sig2.s_dom_blocks
 Sig.s_dom_blocks
 Sig.s_plot


%% check for unvalid Sigma on diagonal

%{
for k = 1:numel(Sigg)
    if ~isempty(Sigg{k})
        negative_diagonal(k) = sum(Sigg{k}.f_diagonal< 0);
        imaginary_diagonal(k) = sum(abs(imag(Sigg{k}.f_diagonal)) > 0);
        disp(['negative diagonal entries: ', num2str(negative_diagonal(k)), ...
        ', imaginary diagonal entries: ', num2str(imaginary_diagonal(k))])
    end
end
%}


%% check circular current

%Sigg{12}.s_plot('A', {min(En.Energies)-10^10, min(En.Energies)+10})

Sig2 = Sigg{22}
Sig2.Values = cellfun(@(x) diag(diag(x)), Sig2.Values, 'uniformoutput', false);
Sig2.Values
A = Sig2.f_full_current(En)

sum(sum(A-A', 2))
% 
%Sigg{12}.s_plot('A', {min(En.Energies)-10^10, min(En.Energies)+10})
    
%}