% TODO optimize the code

function plot_hexagon(a1, a2, line_color, scale, mirrored, gamma_center)
% PLOT_HEXAGON Summary of this function goes here
% plot graphene like hexagon patterns
% Detailed explanation goes here
% a1, a2 unit vector of primitive cells
% scale ~number of hexagons to plot
% mirrored to mirror the pattern against origin
% gamma center if true, the origin is gamma point of the lattice
% this is useful when plotting reciprocal space
% the origin is the center of a hexagon
    grid_start = -10 * scale;
    grid_end = 10 * scale;
    ab_line1 = 1/3 * (a1 + a2);
    ab_line2 = ab_line1 - a1;
    ab_line3 = ab_line1 - a2;
    function plot_unit_cell(center)
        c1 = mirrored * [center, center + ab_line1];
        c2 = mirrored * [center, center + ab_line2];
        c3 = mirrored * [center, center + ab_line3];
        line(c1(1, :), c1(2, :), 'Color', line_color);
        line(c2(1, :), c2(2, :), 'Color', line_color);
        line(c3(1, :), c3(2, :), 'Color', line_color);
    end
    [x, y] = meshgrid(grid_start:grid_end);
    xy = [x(:), y(:)];
    for c_ind = 1:size(xy, 1)
        center = [a1, a2] * xy(c_ind, :)';
        if gamma_center
            center = center + ab_line1;
        end
        plot_unit_cell(center)
    end
    axis equal
end

