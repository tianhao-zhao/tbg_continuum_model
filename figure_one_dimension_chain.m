k = -4*pi:0.01:4*pi;
unit_cell_size = 3;
E = one_dimension_chain(unit_cell_size, k);
fig = figure();
ax = axes(fig);

for i = 1:unit_cell_size
    plot(ax, k, E(i, :));
    hold(ax, 'on');
end
axis(ax, 'equal');

