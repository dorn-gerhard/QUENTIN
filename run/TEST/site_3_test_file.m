%% Setup part 
GROUNDPATH = '/temp/dorn/run/bm_program/';

	addpath(GROUNDPATH)
	addpath([GROUNDPATH, 'classes/'])
	addpath([GROUNDPATH, 'functions/'])
	addpath([GROUNDPATH, 'functions/class_functions/'])
    addpath([GROUNDPATH, 'functions/c_functions/m_GF/'])
	addpath([GROUNDPATH, 'functions/c_functions/'])
	addpath([GROUNDPATH, 'functions/plot/'])
    addpath([GROUNDPATH, 'functions/test_functions/'])
    
    
% run setup with the following settings (energy cut, delta_E, ...)

setup.research_question = '3 site, changed onsite-energy - test sec approx and hope it to fail';
ID = setup_file_site_3(setup);
setup_file = ['./data/setup_', ID, '.mat'];

%ID = 'IRLM_180929_0';


%% Execution part
ID = 'site_3_181206_9';

load(setup_file)

for k = 1:31
    exec_bm(setup_file, num2str(k-1));
    
    %load(['./data/data_', ID, '_', num2str(k-1), '.mat'])
   
    %data.status_full
    %analyse convergence criteria, K ,etc...
    
    
end 
    

%% postprocess

post_run_fuse_files(ID)
%% check setup_file
load(setup_file)


%% analyse


date = '1908';
listing = dir(['./data/combined_data_site_3', '_', date, '*.mat']);
for k = 1:numel(listing)
    disp(listing(k).name)
end
file_name_result = ['./data/', listing(end).name];
load(file_name_result)

%data_combined = transform_data_combined(data_combined, 'diagonal');

print_now = false;



struct2vars(data_combined)
    struct2vars(parameters)


%

II(abs(II) > 10) = nan;
NN(abs(NN) > 30) = nan;

h = 6.626*10^-34;
e = 1.602*10^-19;
e_factor =  2*e^2/h*10^6;

I = II*e_factor;
I_cpt = II_cpt * e_factor;
I_bmcpt = II_bmcpt * e_factor;


num_Vb = numel(Vb_vec);



VB = repmat(Vb_vec, num_Vg,1);
VG = repmat(Vg_vec.', 1, num_Vb);



Vb_delta = VB(2)- VB(1);
VB_delta = ones(size(VB))*Vb_delta;

    
figure
subplot(3,1,1)
plot_settings

pc = plot(VB,I(:,:,2),'-x');
hold on
pc = plot(VB,I(:,:,2)','-x');
pc = plot(VB,I(:,:,3),'-x');
pc = plot(VB,I(:,:,4),'-x');

plot_supertitle(['$\verb|', condor.setup_id, '|$'])
title('BM current')
xlabel('V_b'), ylabel('I')
hold on
grid on

subplot(3,1,2)

pc = plot(VB,-I_cpt(:,:,1),'-x');
hold on
pc = plot(VB,I_cpt(:,:,2),'-x');
pc = plot(VB,-I_cpt(:,:,3),'-x');
pc = plot(VB,I_cpt(:,:,4),'-x');


title('CPT+BMCPT Current')
xlabel('V_b'), ylabel('I')
grid on



pc = plot(VB,-I_bmcpt(:,:,1),'--s');
hold on
pc = plot(VB,I_bmcpt(:,:,2),'--s');
pc = plot(VB,-I_bmcpt(:,:,3),'--s');
pc = plot(VB,I_bmcpt(:,:,4),'--s');

subplot(3,1,3)
xweight = cellfun(@(x) sum(x.Tab.XWeight), Sigg);
pc = plot(VB, xweight, '-.o');
title('Off diagonal Weight')
xlabel('Off diagonal Weight'), ylabel('I')
grid on

if print_now
    plot_pdf(condor.setup_id, parameters)
end


%
% conductivity
%{
dif_I_bmcpt = diff(I_bmcpt(:,:,2), [],2)./Vb_delta;
dif_I_cpt = diff(I_cpt(:,:,2), [],2)./Vb_delta;

dif_I_bm = diff(I(:,:,2), [],2)./Vb_delta;
dif_I_bmcpt = [dif_I_bmcpt(:,1), (dif_I_bmcpt(:,1:end-1) + dif_I_bmcpt(:,2:end))/2];


figure
plot_settings
subplot(3,1,1)

pc = plot(VB(:,1:end-1),dif_I_bm, 'x-');
title([condor.setup_id,' BMConductivity'])
xlabel('V_b'), ylabel('I')
grid on

subplot(3,1,2)

pc = plot(VB(:,1:end-1),dif_I_cpt, '-x');
title('CPTConductivity')
xlabel('V_b'), ylabel('I')
grid on

subplot(3,1,3)

pc = plot(VB(:,1:end-1),dif_I_bmcpt, '-x');
title('BMCPTConductivity')
xlabel('V_b'), ylabel('I')
grid on
%plot_pdf([condor.setup_id, '_conductivity'])

%}
% analysis of data quality

figure
plot_settings
subplot(3,1,1)
plot_supertitle(['$\verb|', condor.setup_id, '|$'])

pc = plot(VB,log10(column_sum), '-x');
title(' Columns sum (logarithmic)')
xlabel('V_b'), ylabel('column\_sum (log10)')
grid on


subplot(3,1,2)

pc = plot(VB,positive_eigenvalues, '-x');
title('Positive Eigenvalues')
xlabel('V_b'), ylabel('Positive Eigenvalues')
grid on

subplot(3,1,3)

pc = plot(VB,multiple_solutions, '-x');
title('Multiple Solutions')
plot_colorbar_integer(8);
xlabel('V_b'), ylabel('Multiple solutions')
grid on
ylim([0,max(multiple_solutions)])
if print_now
    plot_pdf(condor.setup_id)
end
% second evaluation


figure
plot_settings
subplot(3,1,1)

pc = plot(VB,NN, '-x');
plot_supertitle(['$\verb|', condor.setup_id, '|$'])
title('Particle sector')
plot_colorbar_integer;
%caxis([-1, 1]*30)
xlabel('V_b'), ylabel('N')
grid on


subplot(3,1,2)

pc = plot(VB,log10(abs(minimal_value)), '-x');
title('minival eigenvalue (log10)')

plot_colorbar;
xlabel('V_b'), ylabel('minival eigenvalue')
grid on

subplot(3,1,3)

pc = plot(VB,log10(abs(eigenvalue_tolerance)), '-x');
title('eigenvalue_tolerance (log10)')
plot_colorbar;
xlabel('V_b'), ylabel('eigenvalue_tolerance')
grid on
if print_now
    plot_pdf(condor.setup_id)
end



% third evaluation

figure
plot_settings
plot_supertitle(['$\verb|', condor.setup_id, '|$'])
plot(VB, data_combined.energy_cut, 'r')
hold on
%legend('maximal_calc_energy')
grid on

xlabel('VB')
ylabel('occupied eigenstates')
title(' occupation of eigenstates in steady state dependend on applied voltage')
%x = colormap(parula(100));

x = colormap(hsv(100))  ;

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
            if isfile(energy.file{1})
                load(energy.file{1})
            else
                En = get_En(setup_file, Vg_vec(1));
            end
            %load(['./data/energy_Vg_', num2str(Vg_vec(1)), '.mat'])
        end
        energies = En.Energies('E', {-inf, data_combined.energy_cut(k)});
        
        for j = 1:numel(energies)
            plot(VB(k), energies(j), 's', 'markersize', 5, 'color', [1,1,1]*0.70)
        end
        for j = 1:length(T.Index)
            
            WW = T.weight(j);
            if T.weight(j) > 1
                warning ('probability > 1!')
            else
                plot(VB(k), EE(j), 'o', 'markersize', 1+7*WW, 'markerfacecolor', x(ceil(WW*length(x)),:), 'color', x(ceil(WW*length(x)),:))
            end
            hold on
        end
        
    end
end
      
colorbar;

if print_now 
    plot_pdf(condor.setup_id)
end


%Sigg{12}.s_plot('A', {min(En.Energies)-10^10, min(En.Energies)+10})

Sig2 = Sigg{22}
Sig2.Values = cellfun(@(x) diag(diag(x)), Sig2.Values, 'uniformoutput', false);
Sig2.Values
A = Sig2.f_full_current(En)

sum(sum(A-A', 2))
% 
%Sigg{12}.s_plot('A', {min(En.Energies)-10^10, min(En.Energies)+10})
    
% 
% En.s_spectrum([],[],1)
% plot_settings
% 
% if print_now 
%     plot_pdf(condor.setup_id)
% end

% Define Research Answer to Research Question
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

