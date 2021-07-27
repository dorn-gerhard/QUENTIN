function rescale_complex_plot(source, eventdata)
ht = guidata(gcf);

fig = gcf;
k_max = numel(fig.Children);

if ht.button.Value
    for k = 1:k_max - 1
        fig.Children(k_max).Children(k).Children(2).Children(end).CData = ht.C_Data_equal{k};
    end
    xy_max = max(ht.x_max,ht.y_max);
    fig.Children(k_max).SelectedTab.Children(3).XTickLabel = {num2str(-xy_max,3), '0', num2str(xy_max,3)};
    fig.Children(k_max).SelectedTab.Children(3).YTickLabel = {num2str(-xy_max,3), '0', num2str(xy_max,3)};
else
    for k = 1:k_max-1
        fig.Children(k_max).Children(k).Children(2).Children(end).CData = ht.C_Data{k};
    end
    fig.Children(k_max).SelectedTab.Children(3).XTickLabel = {num2str(-ht.x_max,3), '0', num2str(ht.x_max,3)};
    fig.Children(k_max).SelectedTab.Children(3).YTickLabel = {num2str(-ht.y_max,3), '0', num2str(ht.y_max,3)};
    
end

    guidata(gcf,ht);