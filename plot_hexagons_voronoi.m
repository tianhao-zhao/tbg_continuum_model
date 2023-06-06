% plot hexagon with voronoi method
% TODO remove extra outer lines 
function plot_hexagons_voronoi(b1_input, b2_input, color, scale)
    b1 = b1_input(:);
    b2 = b2_input(:);
    grid_start = scale * (-10);
    grid_end = scale * (10);
    [x, y] = meshgrid(grid_start:grid_end);
    xy = [x(:), y(:)]';
    xyt = [b1, b2] * xy;
    xt = xyt(1, :);
    yt = xyt(2, :);
    [xv, yv] = voronoi(yt(:), xt(:));
    ax = axes(figure());
    plot(ax, xv, yv, '-', Color=color);
    axis equal;
    xmin = grid_start * norm(b1);
    xmax = -xmin;
    ymin = grid_start * norm(b1);
    ymax = -ymin;
    axis(ax, [xmin, xmax, ymin, ymax]);
end

