function plot_single_hexagon(ax, center, l, r, highlight, linewidth)
% plot a hexagon in the ceter with l edge length
% rotated by r (radius degree)
    center = center(:);
    rm = [[cos(r), -sin(r)]; [sin(r), cos(r)]];
    a60 = pi / 3;
    r60 = [[cos(a60), -sin(a60)]; [sin(a60), cos(a60)]];
    % line start/end relative to center
    ls = rm*[l; 0];
    le = rm*[l/2; sqrt(3)*l/2];
    lp = [ls, le];
    for ii = 1:6
        line(ax, lp(1, :) + center(1), lp(2, :) + center(2), 'Color', 'black',...
            'LineWidth', linewidth); 
        hold(ax, 'on');
        lp = r60 * lp;
    end
    if highlight
        for ii = 1:6
            if mod(ii, 2) == 1
                scatter(ax, lp(1, 1) + center(1), lp(2, 1) + center(2), 72*linewidth,...
                    'MarkerFaceColor', [1, 0, 0], 'MarkerEdgeColor', 'none');
                hold(ax, 'on');
                scatter(ax, lp(1, 2) + center(1), lp(2, 2) + center(2), 72*linewidth,...
                    'MarkerFaceColor', [0.5, 0, 0], 'MarkerEdgeColor', 'none');
                hold(ax, 'on');
            end
            lp = r60 * lp;
        end
    end

    axis(ax, 'equal');
end

