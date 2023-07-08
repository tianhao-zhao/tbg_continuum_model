k = 0:0.01:2*pi;
for unit_cell_size = 1:10
    E = one_dimension_chain(unit_cell_size, k);
    fig = figure();
    ax = axes(fig);
    for i = 1:unit_cell_size
        plot(ax, k, E(i, :));
        hold(ax, 'on');
    end
    axis(ax, 'equal');
    xline(ax, 0:(2*pi)/unit_cell_size:2*pi, '-');
end