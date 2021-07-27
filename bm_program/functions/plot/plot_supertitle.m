function [] = plot_supertitle(notation, fig)
if nargin < 2 || isempty(fig), fig = gcf; end


plot_settings()

pos = get(gcf,'position'); 
sizetitle = 15; sizelabel = 13;
annotation(fig, 'textbox',  'units','pixel','position', [10, pos(4)*0.035, pos(3)-20,30], 'string', notation, 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center','fontsize', sizelabel , 'interpreter', 'latex')
