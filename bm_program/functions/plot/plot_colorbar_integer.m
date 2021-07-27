function colmap = plot_colorbar_integer(n_index)
if nargin < 1 || isempty(n_index), n_index = 64; end

[~, cmax] = caxis;
if cmax > 30, error('not defined for higher integers'); end

%n_colors = ceil(abs(cmax - cmin))*2+1;

% for testing
%n_index = 128; 
%cmax = 30;

%limits of the input data (0 bis n_integer)

color_range = [0, 0, 0; ...
               0.5, 0, 0; ...
               0, 0, 0.5; ...
               0,   0.5, 0; ...
               0,   0, 1; ...
               0, 0.5, 0.5; ...
               0 , 1, 0; ...
               0.5, 0.5, 0; ...
               1,   0, 0; ...
               0.5, 0, 0.5; ...
               0.5, 0, 1; ...
               0,   0.5,   1; ...
               0,   1,   0.5; ...
               0.5, 1,   0; ...
               1,   0.5,   0; ...
               1,   0,   0.5; ...
               0.5, 0.5, 0.5];
color_range = flipud([color_range; 1- flipud(color_range(1:end-1,:))]);           

%delete some colors you don't like
color_range([1,3:4],:) = [];


integer_max = ceil(cmax);


colors = color_range(1:integer_max+1, :);

delta = ceil(n_index/integer_max);
n_index = delta * integer_max;

colmap = zeros(n_index,3);
for k = 1:3
    for l = 0:integer_max-1
        colmap((1:delta)+delta*l,k) = linspace(colors(l+1,k), colors(l+2,k),delta);
        
    end
end

% for testing
% x = linspace(0,cmax,n_index);
% [xx,yy] = meshgrid(x)
% pcolor(xx,yy,xx)

colorbar
colormap(gca,colmap)

