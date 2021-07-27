function [] = plot_spy(K, sensitivity)
if nargin < 2 || isempty(sensitivity), sensitivity = 10.^(-12:-3); end



for k = 1:numel(sensitivity)-1
    spy(abs(K) >= sensitivity(k) & abs(K) < sensitivity(k+1))
    hold on
end
spy(abs(K) > sensitivity(end))

colmap = parula(numel(sensitivity));


x=get(gca,'children');
for k = 1:numel(sensitivity)
    set(x(end-k+1),'color',colmap(k,:))
end

colorbar
caxis([min(sensitivity), max(sensitivity)])
set(gca,'colorscale','log')