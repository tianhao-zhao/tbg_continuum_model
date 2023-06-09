function add_angle_indicator(ax, center, angle_start, angle_end, radius,...
    str, color, text_shift)
% this function plots an arc as an angle indicator
    center = center(:);
    step = 1/500 * pi;
    if angle_start < angle_end
        angles = angle_start:step:angle_end;
    else
        angles = angle_end:step:angle_start;
    end

    x = radius * cos(angles) + center(1);
    y = radius * sin(angles) + center(2);
    plot(ax, x, y, 'Color', color);
    angle_mid = (angle_start + angle_end)/2;
    text(ax, radius*cos(angle_mid) + center(1) + text_shift(1),...
        radius*sin(angle_mid) + center(2) + text_shift(2), str, 'Color', color);
    hold(ax, 'on');
end

