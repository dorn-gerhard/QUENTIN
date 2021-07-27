function colmap = plot_colorbar2(n_index, neg_color, pos_color, nobar)
if nargin < 1 || isempty(n_index), n_index = 64; end
if nargin < 2 || isempty(neg_color), neg_color = [255, 255, 0]/255; end
if nargin < 3 || isempty(pos_color), pos_color = [125, 46, 143]/255; end
if nargin < 4 || isempty(nobar), nobar = false; end

%colormap - setup
zero_color = [1, 1, 1];

% TODO: add intermediate steps (at first and fourth quantil)
% find colors
third_color = [239, 0, 241]/255;
first_color = [255, 209, 84]/255;


[cmin, cmax] = caxis;
range = linspace(cmin,cmax,n_index);
index_zero = find(min(abs(range)) == abs(range),1)-1;
first_quantil = floor(index_zero/2);
third_quantil = index_zero + floor((n_index - index_zero)/2);

if index_zero < n_index/6 || index_zero > n_index*5/6, 
    index_zero = floor(n_index/2);
    first_quantil = floor(n_index/4);
    third_quantil = floor(n_index/4*3);
end

bars = [first_quantil, index_zero - first_quantil, third_quantil - index_zero, n_index - third_quantil];
colmap = zeros(n_index,3);
for k=1:3
    colmap(:,k) = [linspace(neg_color(k), first_color(k),bars(1)), ...
        linspace(first_color(k), zero_color(k),bars(2)), ...
        linspace(zero_color(k), third_color(k),bars(3)), ...
        linspace(third_color(k), pos_color(k),bars(4))]';
end

if ~nobar
    colorbar
    colormap(gca,colmap)
end

%colormapeditor
