k = 0:0.01:2*pi;
for unit_cell_size = 1:5
    o = one_dimension_chain(unit_cell_size, k);
    E = o.E;
    fig = figure();
    ax = axes(fig);
    for i = 1:unit_cell_size
        plot(ax, k, E(i, :));
        hold(ax, 'on');
    end
    axis(ax, 'equal');
    xline(ax, 0:(2*pi)/unit_cell_size:2*pi, '-');
    V = o.V;
    fig2 = figure();
    ax2 = axes(fig2);
    for i = 1:unit_cell_size
        plot(ax2, V(:, i), 'o');
        hold(ax2, 'on');
    end
end