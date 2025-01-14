% plot hexagon with voronoi method
% it is aimed to replace plot_hexagons 
% will add customized shift first, then add gamma shift, then mirror the
% lattice
function plot_hexagons_voronoi(ax, a1, a2, line_color, scale, mirrored,...
    gamma_center, nvargs)
    arguments
        ax;
        a1;
        a2;
        line_color;
        scale (1, 1) {mustBeInteger};
        mirrored;
        gamma_center;
        nvargs.shift (2, 1) double = [0; 0];
        nvargs.highlight (1, 1) logical = 0;
        nvargs.lw (1, 1) double = 1.5
    end
    if line_color == 'r'
        line_color = [1, 0, 0];
    end
    if line_color == 'b'
        line_color = [0, 0, 1];
    end
    
    unit_cell_cnt = 2;
    a1 = a1(:);
    a2 = a2(:);
    grid_start = scale * (-unit_cell_cnt);
    grid_end = scale * (unit_cell_cnt);
    [x, y] = meshgrid(grid_start:grid_end);
    xy = [x(:), y(:)]';
    xyt = [a1, a2] * xy;
    xt = xyt(1, :) + nvargs.shift(1);
    yt = xyt(2, :) + nvargs.shift(2);
    if ~gamma_center
        gamma_shift = 1/3 * (a1 + a2);
        xt = xt - gamma_shift(1);
        yt = yt - gamma_shift(2);
    end
    xt = mirrored * xt;
    yt = mirrored * yt;
    % try using voronoi limit
    % bs_ext = unit_cell_cnt*[-a1-a2, a1-a2, a1+a2, -a1+a2]';
    % [V,C,XY]=VoronoiLimit(xt,yt,'bs_ext',bs_ext);
    % xv = V(1, :);
    % yv = V(2, :);
    [xv, yv] = voronoi(xt(:), yt(:));
    plot(ax, xv, yv, '-', Color=line_color, LineWidth=nvargs.lw);
    hold(ax, 'on');
    A_atoms = [a1, a2] * xy + nvargs.shift;
    ab_line = 1/3 * (a1 + a2);
    B_atoms = [a1, a2] * xy + nvargs.shift + ab_line;
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
    xmin = grid_start * norm(a1);
    xmax = -xmin;
    ymin = grid_start * norm(a2);
    ymax = -ymin;
    axis(ax, [xmin, xmax, ymin, ymax]);
    hold(ax, 'on');
end

