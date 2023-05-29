% plot raw lattice with all unit vectors connected
function plot_raw_lattice(b1_input, b2_input, color, scale)
    b1 = transpose(b1_input);
    b2 = transpose(b2_input);
    grid_start = scale * (-10);
    grid_end = scale * (10);
    [x, y] = meshgrid(grid_start:grid_end);
    xy = [x(:), y(:)];
    T = [b1;b2];
    xyt = xy * T;
    xt = reshape(xyt(:,1),size(x));
    yt = reshape(xyt(:,2),size(y));
    plot(xt,yt,'ro-', color=color);
    hold on;
    plot(xt',yt','r-', color=color);
    axis equal;
    %axis square
end

