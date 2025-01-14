% TODO optimize the code
% see plot_hexagon_voronoi for a second method
function plot_hexagons(ax, a1, a2, line_color, scale, mirrored,...
    gamma_center, nvargs)
% PLOT_HEXAGON Summary of this function goes here
% plot graphene like hexagon patterns
% Detailed explanation goes here
% a1, a2 unit vector of primitive cells
% scale ~number of hexagons to plot
% mirrored to mirror the pattern against origin
% gamma center if true, the origin is gamma point of the lattice
% this is useful when plotting reciprocal space
% the origin is the center of a hexagon
% highlight if true, plot dot-like atoms
    arguments
        ax;
        a1;
        a2;
        line_color;
        scale;
        mirrored;
        gamma_center;
        nvargs.shift (1, 2) double = [0.0, 0.0];
        nvargs.highlight (1, 1) logical = 0;
        nvargs.lw (1, 1) double = 1.5;
    end

    % temporary code, convert color string into vector
    if line_color == 'r'
        line_color = [1, 0, 0];
    end
    if line_color == 'b'
        line_color = [0, 0, 1];
    end
      
      
    nvargs.shift = nvargs.shift(:);
    grid_start = -2 * scale;
    grid_end = 2 * scale;
    ab_line1 = 1/3 * (a1 + a2);
    ab_line2 = ab_line1 - a1;
    ab_line3 = ab_line1 - a2;
    function plot_unit_cell(ax, center)
        c1 = mirrored * [center, center + ab_line1];
        c2 = mirrored * [center, center + ab_line2];
        c3 = mirrored * [center, center + ab_line3];
        line(ax, c1(1, :), c1(2, :), 'Color', line_color, LineWidth=nvargs.lw);
        line(ax, c2(1, :), c2(2, :), 'Color', line_color, LineWidth=nvargs.lw);
        line(ax, c3(1, :), c3(2, :), 'Color', line_color, LineWidth=nvargs.lw);
        hold(ax, 'on');
    end
    [x, y] = meshgrid(grid_start:grid_end);
    xy = [x(:), y(:)];
    for c_ind = 1:size(xy, 1)
        center = [a1, a2] * xy(c_ind, :)' + nvargs.shift;
        if gamma_center
            center = center + ab_line1;
        end
        plot_unit_cell(ax, center);
    end
    hold(ax, 'on')
    A_atoms = [a1, a2] * xy' + nvargs.shift;
    B_atoms = A_atoms + ab_line1(:);
    if nvargs.highlight
        scatter(ax, A_atoms(1, :), A_atoms(2, :), 72*nvargs.lw,...
            'MarkerFaceColor',line_color,...
            'MarkerEdgeColor','none');
        hold(ax, 'on');
        scatter(ax, B_atoms(1, :), B_atoms(2, :), 72*nvargs.lw,...
            'MarkerFaceColor',[0.7, 0.7, 0.7]+0.3*line_color,...
            'MarkerEdgeColor','none');
    end
    axis(ax, 'equal');
    hold(ax, 'on');
end

