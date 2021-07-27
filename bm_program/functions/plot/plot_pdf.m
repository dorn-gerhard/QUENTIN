function plot_pdf(file_name, parameters, crop_true)
%TODO: define a standardized way to add properties to standard properties
if nargin < 3 || isempty(crop_true), crop_true = true; end
if nargin < 2 || isempty(parameters), parameters = []; end
crop_tag = [];
if crop_true
    crop_tag = '-nocrop'; 
end


export_fig([file_name, '.pdf'], '-pdf', '-q101', '-painters', '-append', crop_tag)

if ~isempty(parameters)
    for k = 1:numel(parameters)
        
        figure
        plot_settings
        text_param = plot_parameters(parameters(k));
        annotation('textbox',[0.05 0.05 1 0.90], 'String', text_param,'FitBoxToText','on', 'fontsize', 8);
        
        export_fig([file_name, '.pdf'], '-pdf', '-q101', '-painters', crop_tag, '-append')
    end
end


end