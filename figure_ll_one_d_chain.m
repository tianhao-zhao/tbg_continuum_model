% This script is to plot the magnetic one dimensional chain

% physics variables
a = 0.01;
Delta_n = 5;
Delta = Delta_n * a;
q = 3;
p = 2;
N = 30;

% plotting variables
ball_size = 1;
text_size = 1;
indication_size = 1;

fig = figure();
ax = axes(fig);

% plot a chain of balls
ball_pos_x = 0:a:(N - 1)*a;
ball_pos_y = zeros(1, N);
color_k = 2*pi*p/(q*Delta);
ball_color_params = sin(color_k*ball_pos_x);
scatter(ball_pos_x, ball_pos_y, [],...
    ball_color_params, 'o', 'filled');
colormap(ax, 'cool');
hold(ax, 'on');

% plot arcs
for i = 1:Delta_n:(N-Delta_n)
    add_arc(ax, ball_pos_x(i), ball_pos_y(i), ...
        ball_pos_x(i+Delta_n), ball_pos_y(i+Delta_n), indication_size, 'black');
end
for i = 2:Delta_n:(N-Delta_n)
    add_arc(ax, ball_pos_x(i), ball_pos_y(i), ...
        ball_pos_x(i+Delta_n), ball_pos_y(i+Delta_n), indication_size, [0.5, 0.5, 0.5]);
end

% add double arrow
annotation(fig, 'doublearrow', [0.03, 0.53], [.46, .46]);
% add text
text(ax, Delta_n*a/2, 1.5*a, '\Delta', FontSize=24*text_size);
text(ax, q * Delta_n * a/2, -1.5*a, sprintf('q = %d, p = %d', [q, p]),...
    FontSize=24*text_size);

% post process
axis(ax, 'off');
axis(ax, 'equal');

tighten_margin(fig, ax)
saveas(fig, '..\..\figures\magnetic_chain.png');