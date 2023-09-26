function add_arc(ax, x1, y1, x2, y2, lw, color, theta)
    arguments
        ax;
        x1;
        y1; 
        x2; 
        y2;
        lw = 1;
        color = 'b';
        theta = pi / 3;
    end
    dist = norm([x2 - x1, y2 - y1]);
    if x2 ~= x1
        tilt_angle = atan((y2 - y1) / (x2 - x1));
    else
        tilt_angle = pi / 2;
    end
    r = (dist/2) / sin(theta/2);
    x_mid = (x1 + x2)/2;
    y_mid = (y1 + y2)/2;
    center_x = x_mid + r*cos(theta/2)*sin(tilt_angle);
    center_y = y_mid - r*cos(theta/2)*cos(tilt_angle);
    if y2 ~= y1
        center_angle = atan((x1 - x2) / (y2 - y1));
    else
        center_angle = pi / 2;
    end

    % number of points
    N = 1000;
    delta_angle = theta / 200;
    angle_list = (center_angle - theta/2):delta_angle:(center_angle + theta/2);
    x_list = center_x + r*cos(angle_list);
    y_list = center_y + r*sin(angle_list);

    % plot
    plot(ax, x_list, y_list, Color=color, LineWidth=lw*0.5);
    hold(ax, 'on');
end