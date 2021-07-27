function   Hub_spectrum(basis, system, E_cut)

%% global parameters
rnd = 10;

% GROUNDPATH = '/temp/dorn/run/bm_program/';
% 
% 	addpath(GROUNDPATH)
% 	addpath([GROUNDPATH, 'classes/'])
% 	addpath([GROUNDPATH, 'functions/'])
% 	addpath([GROUNDPATH, 'functions/class_functions/'])
% 	addpath([GROUNDPATH, 'functions/c_functions/'])
% 	addpath([GROUNDPATH, 'functions/plot/'])


h.fig = figure('position',[100,100,1400,900]);
set(h.fig, 'color','w')
% h.ax1 = subplot(3,3,1);
% line(0,0)
% h.ax2 = subplot(3,3,2);
% line(0,0)
% h.ax3 = subplot(3,3,3);
% line(0,0)
% h.ax4 = subplot(3,3,4);
% line(0,0)
% h.ax5 = subplot(3,3,5);
% line(0,0)
 h.ax6 = subplot(2,2,1);

 
 
 
 
 
%% variables create
if nargin < 3 || isempty(E_cut)
    h.E_cut = inf;
else
    h.E_cut = E_cut;
end

if nargin < 2 || isempty(system)
    h.U = 0;
    h.V = 0;
    h.b = 0;
    h.N_site = 3;
    h.system.xi =0;
    h.system.U = 0;
    h.system.V = 0;
    h.system.b = 1;
    h.system.part_hole_sym = true;

else
    h.system.xi =system.xi;
    h.system.U = system.U;
    h.system.V = system.V;
    h.system.b = system.b;
    h.system.part_hole_sym = system.part_hole_sym;

    h.U = system.U;
    h.V = system.V; 
    h.b = system.b;
    h.N_site = length(system.U);
end
    
h.xi0 = 0;
h.xi = 0;
h.xi3 = 0;
h.chain = false;
if nargin < 1 || isempty(basis)
    h.basis.N_sector = [];
    h.basis.S_sector = [];
    h.basis.Docc_input = 1;
    h.basis.N_spin = 2;
    h.basis.N_site = 2;
else
    h.basis.N_sector = basis.N_sector;
    h.basis.S_sector = basis.S_sector;
    h.basis.Docc_input = basis.Docc_input;
    h.basis.N_spin = basis.N_spin;
    h.basis.N_site = basis.N_site;
end



h.geometry.structure = 'p';
h.Geom = Geometry(h.basis.N_site, h.geometry.structure, ones(h.basis.N_site));
h.Bas = Basis(h.basis.N_site, h.basis.N_sector, h.basis.S_sector, h.basis.Docc_input, h.basis.N_spin); % class
        

h.energy.calc.tolerance = 10^-14;
h.energy.calc.max_number = 10;
h.energy.second_cut = 10;

guidata(h.fig, h)

%% uicontrols create
h = guidata(gcf);
sl_len = 300;
sl_hei = 30;
text_len = 80;
gap = 7;

h.slider_xi = uicontrol('style', 'slider', 'min',-5, 'max', 5, 'value',0, 'position', [10,30,sl_len,sl_hei], 'sliderstep', [0.01,0.1]);
h.sl_lab_xi = uicontrol('style', 'text', 'position', [sl_len + 60, 30-gap, text_len, sl_hei], 'string', ['xi: ', num2str(0)], 'backgroundcolor','w');
h.slider_U = uicontrol('style', 'slider', 'min',-5, 'max', 10, 'value',0, 'position', [10,70,sl_len,sl_hei], 'sliderstep', [0.01,0.1]);
h.sl_lab_U = uicontrol('style', 'text', 'position', [sl_len + 60, 70-gap, text_len, sl_hei], 'string', ['U: ', num2str(0)], 'backgroundcolor','w');


h.slider_b = uicontrol('style', 'slider', 'min',-5, 'max', 5, 'value',0, 'position', [10,110,sl_len,sl_hei], 'sliderstep', [0.01,0.1]);
h.sl_lab_b = uicontrol('style', 'text', 'position', [sl_len + 60, 110-gap, text_len, sl_hei], 'string', ['b: ', num2str(0)], 'backgroundcolor','w');
h.slider_V = uicontrol('style', 'slider', 'min',-5, 'max', 6, 'value',0, 'position', [10,150,sl_len,sl_hei], 'sliderstep', [0.01,0.1]);
h.sl_lab_V = uicontrol('style', 'text', 'position', [sl_len + 60, 150-gap, text_len, sl_hei], 'string', ['V: ', num2str(0)], 'backgroundcolor','w');

h.sl_lab_xi3 = uicontrol('style', 'text', 'position', [sl_len + 120, 30-gap, text_len, sl_hei], 'string', ['xi_3: ', num2str(0)], 'backgroundcolor','w');
h.slider_xi3 = uicontrol('style', 'slider', 'min',-2, 'max', 2, 'value',0, 'position', [sl_len + 200, 30 - gap ,text_len,sl_hei], 'enable', 'inactive', 'sliderstep', [0.02,0.2]);

h.button_chain = uicontrol('style','togglebutton',...
                        'position', [2.3*sl_len, 30, text_len, sl_hei], ...
                        'string', 'Chain');
                    
h.dropdown_N_site = uicontrol('style', 'popupmenu', 'position', [2.3*sl_len, 90, text_len, sl_hei], ...
                        'string', {'1','2','3','4','5','6'});
if nargin >= 1 && ~isempty(basis)
    h.dropdown_N_site.Value = basis.N_site;
end
                    
h.dd_lab = uicontrol('style', 'text', 'position', [2.3*sl_len , 125, text_len, 20], 'string', 'N_site', 'backgroundcolor','w');

h.table_U = uitable('Parent', h.fig,...
                    'Data', h.system.U, ...
                    'Position', [20,300,300,160],...
                    'ColumnWidth', {40,40,40,40,40,40},...
                    'ColumnEditable', true, ...
                    'tooltip', sprintf('U \n set the values \n for the interaction between \n electrons of different spins'),...
                    'CellEditCallback', @set_U_callback);
                
h.table_b = uitable('Parent', h.fig,...
                    'Data', h.system.b, ...
                    'Position', [400,300,300,160],...
                    'ColumnWidth', {40,40,40,40,40,40},...
                    'ColumnEditable', true, ...
                    'tooltip', sprintf('b \n set the values \n for the hopping and \n on-site energies \nTODO: for different spins'),...
                    'CellEditCallback', @set_b_callback);
                
                
h.table_V = uitable('Parent', h.fig,...
                    'Data',h.system.V, ...
                    'Position', [800,300,300,160],...
                    'ColumnWidth', {40,40,40,40,40,40},...
                    'ColumnEditable', true, ...
                    'tooltip', sprintf('b \n set the values \n for the hopping and \n on-site energies \nTODO: for different spins'),...
                    'CellEditCallback', @set_V_callback);
                
h.button_autoscale = uicontrol('style','togglebutton',...
                        'position', [30, 500, text_len, sl_hei], ...
                        'string', 'Autoscale', ...
                        'callback', @set_autoscale_callback);
                
h.button_spinless = uicontrol('style','togglebutton',...
                        'position', [30, 550, text_len, sl_hei], ...
                        'string', 'Spinless', ...
                        'callback', @set_spinless_callback);

h.button_instant = uicontrol('style','togglebutton',...
                        'position', [30, 600, text_len, sl_hei], ...
                        'string', 'Instant calculation', ...
                        'callback', @set_instant_callback);
                    
h.button_calcnow = uicontrol('style','pushbutton',...
                        'position', [30, 650, text_len, sl_hei], ...
                        'string', 'Calculate now', ...
                        'callback', @set_N_site_callback);
                
%uitable(h.fig,'Data',[1,0;0,1]);
%


%% uicontrols initialize

set(h.fig,'resizefcn', {@positioning})
set(h.fig,'WindowButtonDownFcn', @clickfnc);

set(h.button_chain, 'callback', {@set_input_data_callback});

set(h.slider_xi, 'callback', {@set_xi_callback});

set(h.slider_U, 'callback', {@set_U_callback});           

set(h.slider_b, 'callback', {@set_b_callback});   

set(h.slider_V, 'callback', {@set_V_callback});

set(h.dropdown_N_site, 'callback', {@set_N_site_callback});

set(h.slider_xi3, 'callback', {@set_xi3_callback});




guidata(h.fig, h)

if nargin >= 2 && ~isempty(system) && ~isempty(basis)
    set_input_data_callback()
end

%% functions

    function set_xi3_callback(source,eventdata)
        ht = guidata(gcf);
        disp(['xi: ' num2str(ht.xi)])
        
        ht.xi3 = round(rnd*get(ht.slider_xi3, 'value'))/rnd;
        %ht.xi3 = str2double(get(ht.edit_xi3, 'string'));
         set(ht.sl_lab_xi3, 'string', ['xi3: ', num2str(ht.xi3)]);
        guidata(gcf,ht);
        set_input_data_callback();
    end


    function set_N_site_callback(source, eventdata)
        ht = guidata(gcf);
        N_site = get(ht.dropdown_N_site, 'value');
        ht.basis.N_site = N_site;
        
        %update Basis and Geometry
        ht.Bas = Basis(ht.basis.N_site, ht.basis.N_sector, ht.basis.S_sector, ht.basis.Docc_input, ht.basis.N_spin); % class
        ht.Geom = Geometry(ht.basis.N_site, ht.geometry.structure, ones(ht.basis.N_site));
        %update input tables
        old_N_site =  length(ht.table_U.Data);
        diff_N_site = N_site - old_N_site;
        if diff_N_site <= 0
            ht.table_U.Data = ht.table_U.Data(1:N_site, 1:N_site);
            ht.table_b.Data = ht.table_b.Data(1:N_site, 1:N_site);
            ht.table_V.Data = ht.table_V.Data(1:N_site, 1:N_site);
        else
            temp_U = ht.table_U.Data;
            temp_U(1:end-diff_N_site, 1:end-diff_N_site) = 0;
            ht.table_U.Data(end+diff_N_site, end+diff_N_site) = 0;
            ht.table_U.Data(diff_N_site+1:end,diff_N_site+1:end) = ht.table_U.Data(diff_N_site+1:end,diff_N_site+1:end) + temp_U;
            
            temp_b = ht.table_b.Data;
            temp_b(1:end-diff_N_site, 1:end-diff_N_site) = 0;
            ht.table_b.Data(end+diff_N_site, end+diff_N_site) = 0;
            ht.table_b.Data(diff_N_site+1:end,diff_N_site+1:end) = ht.table_b.Data(diff_N_site+1:end,diff_N_site+1:end) + temp_b;
            
            temp_V = ht.table_V.Data;
            temp_V(1:end-diff_N_site, 1:end-diff_N_site) = 0;
            ht.table_V.Data(end+diff_N_site, end+diff_N_site) = 0;
            ht.table_V.Data(diff_N_site+1:end,diff_N_site+1:end) = ht.table_V.Data(diff_N_site+1:end,diff_N_site+1:end) + temp_V;
            
        end
        ht.system.U = ht.table_U.Data;
        ht.system.b = ht.table_b.Data;
        ht.system.V = ht.table_V.Data;
        guidata(gcf, ht);
        set_input_data_callback();
    end

    function set_autoscale_callback(source, eventdata)
        ht = guidata(gcf);
        if ht.button_autoscale.Value
            ht.ax6.XLimMode = 'auto';
            ht.ax6.YLimMode = 'auto';
        else
            ht.ax6.XLimMode = 'manual';
            ht.ax6.YLimMode = 'manual';
        end
    end
        
    function set_spinless_callback(source, eventdata)
        ht = guidata(gcf);
        
        if ht.button_spinless.Value
            ht.basis.N_spin = 1;
        else
            ht.basis.N_spin = 2;
        end
        guidata(gcf,ht);
        set_N_site_callback()
    end

    function set_instant_callback(source, eventdata)
       
    end
        
        
    function set_xi_callback(source, eventdata)
        ht = guidata(gcf);
        ht.xi0 = round(rnd*get(ht.slider_xi, 'value'))/rnd;
        set(ht.sl_lab_xi, 'string', ['xi: ', num2str(ht.xi0)]);
        guidata(gcf, ht);
        set_input_data_callback();
    end

    function set_U_callback(source, eventdata)
        ht = guidata(gcf);
        %change to data from table
        ht.system.U = ht.table_U.Data;
        disp('update_U')
        %ht.system.U = round(rnd*get(ht.slider_U, 'value'))/rnd;
        %set(ht.sl_lab_U, 'string', ['U: ', num2str(ht.U)]);
        guidata(gcf, ht);
        if ht.button_instant.Value
            set_input_data_callback();
        end
    end
    function set_b_callback(source, eventdata)
        ht = guidata(gcf);
        ht.system.b = ht.table_b.Data;
        %ht.system.b = round(rnd*get(ht.slider_b, 'value'))/rnd;
        %set(ht.sl_lab_b, 'string', ['b: ', num2str(ht.b)]);
        guidata(gcf, ht);
        if ht.button_instant.Value
            set_input_data_callback();
        end
    end
    function set_V_callback(source, eventdata)
        ht = guidata(gcf);
        ht.system.V = ht.table_V.Data;
        %ht.system.V = round(rnd*get(ht.slider_V, 'value'))/rnd;
        %set(ht.sl_lab_V, 'string', ['V: ', num2str(ht.V)]);       
        guidata(gcf, ht);
        if ht.button_instant.Value
            set_input_data_callback();
        end
    end

    function clickfnc(source, eventdata)
        
    end

    function set_input_data_callback(source,eventdata)
        ht = guidata(gcf);
        
        ht.chain = get(ht.button_chain,'value');
        if ht.N_site > 8
            set(ht.slider_xi3, 'enable', 'on')
            disp('set to enable')
            ht.system.xi = ones(1,ht.N_site)*ht.xi0;
            ht.system.xi(3) = ht.xi(3) + ht.xi3;
        end
        
        
        
        
        ht.Hub = Hubbard(ht.Geom, ht.system.b , ht.system.U, ht.system.xi, ht.system.V, ht.Bas, ht.system.part_hole_sym); % class

        ht.En = ht.Hub.f_solve(ht.energy.calc, ht.energy.second_cut);
        guidata(gcf, ht);
        xlimit = ht.ax6.XLim;
        ylimit = ht.ax6.YLim;
        details = true;
        ht.En.s_spectrum([],[],details,ht.ax6);
        
        if strcmp( ht.ax6.XLimMode, 'manual')
            xlim(xlimit)
        end
        if strcmp(ht.ax6.YLimMode, 'manual')
            ylim(ylimit)
        end
    end


    function plot_all(source,eventdata)
        ht = guidata(gcf);

        %plot spektrum
        axes(ht.ax6)
        marker = {'x','<','^','s','p','h'};
        

        
        for k = 1:size(ht.ew_u,1)
            if ht.ew_u(k,4) < 7
                marker_plot = marker{ht.ew_u(k,4)};
            else marker_plot = '*';
            end
            plot(ht.ew_u(k,1), ht.ew_u(k,2),marker_plot);
            hold on
        end
        pl = plot(nan,nan, 'bx', nan,nan,'<b', nan, nan,'^b', nan, nan, 'bs', nan, nan,'bp', nan, nan, 'bh', nan,nan,'b*');
        hold off
        legend(pl, '1x','2x','3x', '4x', '5x', '6x', 'more','location','northeastoutside')
        xlabel('N_el')
        ylabel('Energy')
        
        guidata(gcf, ht);
    end

    function positioning(source,eventdata)
        %     ht = guidata(gcf);
        %
        %     window_size = get(ht.figure1,'position');
        %     height = window_size(4);
        %     width = window_size(3);
        %     min_gap = 10;
        %     slider_height = 18;
        %     slider_length = 90;
        %     button_height = 25;
        %     button_length = 80;
        %
        %
        %     set(ht.button1, 'units', 'pixels', 'position', [min_gap, min_gap, button_length, button_height])
        %     set(ht.button2, 'units', 'pixels', 'position', [button_length + 2*min_gap, min_gap, button_length, button_height])
        %     set(ht.toggle_single_patch, 'units', 'pixels', 'position', [2*button_length + 3*min_gap, min_gap, button_length, button_height])
        %
        %     set(ht.slider1, 'units', 'pixels', 'position', [min_gap, height-50, button_length, slider_height])
        %     set(ht.slider2, 'units', 'pixels', 'position', [slider_length + 2*min_gap, height-50, button_length, slider_height])
        %     set(ht.slider3, 'units', 'pixels', 'position', [2*slider_length + 3*min_gap, height-50, button_length, slider_height])
        %     set(ht.slider4, 'units', 'pixels', 'position', [3*slider_length + 4*min_gap, height-50, button_length, slider_height])
        %
        %     set(ht.text1, 'units', 'pixels', 'position', [min_gap, height-70, button_length, slider_height])
        %     set(ht.text2, 'units', 'pixels', 'position', [slider_length + 2*min_gap, height-70, button_length, slider_height])
        %     set(ht.text3, 'units', 'pixels', 'position', [2*slider_length + 3*min_gap, height-70, button_length, slider_height])
        %     set(ht.text4, 'units', 'pixels', 'position', [3*slider_length + 4*min_gap, height-70, button_length, slider_height])
        %
        %
        %     set(ht.text_frequency, 'units', 'pixels', 'position', [min_gap, height-30, button_length, slider_height])
        %     set(ht.text_c_min, 'units', 'pixels', 'position', [slider_length + 2*min_gap, height-30, button_length, slider_height])
        %     set(ht.text_c_max, 'units', 'pixels', 'position', [2*slider_length + 3*min_gap, height-30, button_length, slider_height])
        %     set(ht.text_patch_number, 'units', 'pixels', 'position', [3*slider_length + 4*min_gap, height-30, button_length, slider_height])
        %     disp((button_height + 2*min_gap)/height)
        %     ax1_size = min(0.9 * width, 0.6 * height);
        %     set(ht.axes1, 'units', 'pixel', 'position', [3*min_gap, (button_height + 4*min_gap), ax1_size, ax1_size])
        %     set(ht.axes2, 'units','pixel','position', [3*min_gap, ax1_size + 10 * min_gap, ax1_size, min(0.2*height, height-5*slider_height-ax1_size-10*min_gap)])
    end
end
