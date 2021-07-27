function [] = plot_GF(GF, Grid,  c_site_1, c_site_2)
plot_all = false;
if nargin < 3, plot_all = true; end
if nargin < 4, plot_all = true; end

if isa(Grid,'Grid')
    omega = Grid.Points;
else
    omega = Grid;
end

if ndims(GF) < 3 %NOTE: 
    [sx,sy] = size(GF);
    if sx == 1 || sy == 1
        GF = reshape(GF(:), 1,1,numel(GF));
    else
        if sx > sy
            GF = reshape(GF, sy,sy,sx/sy);
        end
        if sy > sx
            GF = reshape(GF, sx,sx,sy/sx);
        end
        
    end
end
if plot_all
    n_site_1 = size(GF,1);
    n_site_2 = size(GF,2);
    fig = figure;
    plot_settings
    tabgrp = uitabgroup(fig);
    sum_tab = uitab(tabgrp, 'Title', 'overlaided G_{ij}', 'backgroundcolor', 'w');
    trace_tab = uitab(tabgrp, 'Title', 'Tr(G_{ij})', 'backgroundcolor', 'w');
    axes('parent', trace_tab)
    set(gcf, 'color', 'w')
    subplot(2,1,1)
    hold on
    plot(omega, real(diagsum(GF,1,2)))
    grid on
    title('real part')
    xlabel('\omega')
    
    subplot(2,1,2)
    hold on
    plot(omega, imag(diagsum(GF,1,2)) )
    grid on
    title('imaginary part')
    xlabel('\omega')
    
    legend_overlay = cell(1,n_site_1*n_site_2);
    index = 1;
    for c_site_1 = 1:n_site_1
        
        for c_site_2 = 1:n_site_2
            current_tab = uitab(tabgrp, 'Title', ['G_{', num2str(c_site_1), num2str(c_site_2), '}'], 'backgroundcolor', 'w');
            axes('parent', current_tab)
            subplot(2,1,1)
            hold on
            plot(omega, real(squeeze(GF(c_site_1, c_site_2, :))))
            grid on
            title('real part')
            xlabel('\omega')
            
            subplot(2,1,2)
            hold on
            plot(omega, imag(squeeze(GF(c_site_1, c_site_2, :))) )
            grid on
            title('imaginary part')
            xlabel('\omega')
            
            axes('parent', sum_tab)
            subplot(2,1,1)
            hold on
            plot(omega, real(squeeze(GF(c_site_1, c_site_2, :))))
            grid on
            title('real part')
            xlabel('\omega')
            
            subplot(2,1,2)
            hold on
            plot(omega, imag(squeeze(GF(c_site_1, c_site_2, :))) )
            grid on
            title('imaginary part')
            xlabel('\omega')
            legend_overlay{index} = ['G_{', num2str(c_site_1), num2str(c_site_2), '}'];
            index = index+1;
        end
    end
    axes('parent', sum_tab)
    subplot(2,1,1)
    legend(legend_overlay)
    subplot(2,1,2)
    legend(legend_overlay)
else
    
    subplot(2,1,1)
    hold on
    plot(omega, real(squeeze(GF(c_site_1, c_site_2, :))))
    grid on
    title('real part')
    xlabel('\omega')
    
    subplot(2,1,2)
    hold on
    plot(omega, imag(squeeze(GF(c_site_1, c_site_2, :))) )
    grid on
    title('imaginary part')
    xlabel('\omega')
    
end

