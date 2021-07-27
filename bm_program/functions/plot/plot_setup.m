function system_tab = plot_setup(setup_ID, print_now, uitabgr)
if nargin < 2 || isempty(print_now), print_now = false; end
if nargin < 3 || isempty(uitabgr)
    fig = figure;plot_settings
    tabgrp = uitabgroup(fig);
else 
    tabgrp = uitabgr;
end
    



load(['./data/setup_',setup_ID])

struct2vars(parameters)


system_tab = uitab(tabgrp, 'Title', 'System', 'backgroundcolor', 'w');
%axes('parent', system_tab)
tabgrp.SelectedTab = system_tab;
annotation(system_tab, 'textbox','position', [0.3,0.92,0.05,0.03], 'string', 'T_{ij}', 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center' , 'interpreter', 'tex', 'fontname', 'DejaVu Sans Mono')
%annotation(system_tab, 'textbox',[0.17,0.82,0.3,0.1], 'String', num2str(system.b,3),'FitBoxToText','on', 'fontsize', 8);
annotation(system_tab, 'textbox',[0.17,0.82,0.3,0.1], 'String', cell2mat(print_aligned_table(system.b,'  ')),'FitBoxToText','on', 'fontsize', 8, 'fontname', 'DejaVu Sans Mono');

annotation(system_tab, 'textbox','position', [0.3,0.8,0.05,0.03], 'string', 'U_{ij}', 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center' , 'interpreter', 'tex', 'fontname', 'DejaVu Sans Mono')
annotation(system_tab, 'textbox',[0.17,0.7,0.3,0.1], 'String', cell2mat(print_aligned_table(system.U,'  ')),'FitBoxToText','on', 'fontsize', 8, 'fontname', 'DejaVu Sans Mono');

annotation(system_tab, 'textbox','position', [0.3,0.68,0.05,0.03], 'string', 'V_{ij}', 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center' , 'interpreter', 'tex', 'fontname', 'DejaVu Sans Mono')
annotation(system_tab, 'textbox',[0.17,0.58,0.3,0.1], 'String', cell2mat(print_aligned_table(system.V,'  ')),'FitBoxToText','on', 'fontsize', 8, 'fontname', 'DejaVu Sans Mono');

% Hubbard Hamiltonian (interactive (Vb + Vg)

% tables (spin dependent)
% limits in system size?
% number of used values? 0 if 0

% show spectrum

% base facts: N_part, N_spin


%% Bath (+Beta, Vbias range - interactive?)


%% Coupling


%% Condor facts


%% Spectrum analysis


%% energy setup


%% 








end