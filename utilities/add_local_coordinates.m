function add_local_coordinates(ax, center, angle, length, color, lw, ls)
    center = center(:);
    x1 = length/2 * cos(angle);
    y1 = length/2 * sin(angle);
    x2 = -length/2 * sin(angle);
    y2 = length/2 * cos(angle);
    xs1 = center(1) - x1;
    xe1 = center(1) + x1;
    ys1 = center(2) - y1;
    ye1 = center(2) + y1;
    xs2 = center(1) - x2;
    xe2 = center(1) + x2;
    ys2 = center(2) - y2;
    ye2 = center(2) + y2;

    quiver(ax, xs1, ys1, xe1 - xs1, ye1 - ys1,...
        "AutoScale", 'off', 'Color', color, 'LineWidth', lw, 'LineStyle', ls);
    hold(ax, 'on');
    quiver(ax, xs2, ys2, xe2 - xs2, ye2 - ys2,...
        "AutoScale", 'off', 'Color', color, 'LineWidth', lw, 'LineStyle', ls);
    hold(ax, 'on');
end

