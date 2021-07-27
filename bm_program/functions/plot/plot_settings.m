function plot_settings(portrait, monitor)
if nargin < 2 || isempty(monitor), monitor = 1; end
if nargin < 1 || isempty(portrait), portrait = true; end
win_size = get(0,'screensize'); %get(0,'MonitorPositions');
set(gcf,'paperunits', 'centimeters')
set(gcf, 'color', 'w')
if portrait
    set(gcf,'position', [win_size(monitor,3)/2, 0.05*win_size(monitor,3), round(1/sqrt(2)*win_size(monitor,4)*0.8), win_size(monitor,4)*0.8])
    set(gcf,'paperorientation', 'portrait')
else
    
    set(gcf,'position', [win_size(monitor,3)/5,  0.05*win_size(monitor,3), round(sqrt(2)*win_size(monitor,4)*0.8), win_size(monitor,4)*0.8])
    set(gcf,'paperorientation', 'landscape')
end
set(gcf, 'papertype', 'a4')
set(gcf, 'paperpositionmode', 'auto')