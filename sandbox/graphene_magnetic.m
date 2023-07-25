% this code is to test 
% graphene landau level in magnetic field 
% doesn't depend on q or k

N = 400;
a = diag(sqrt(1:N-1), 1);
ad = diag(sqrt(1:N-1), -1);
EA = 10;
EB = 10;

B_max = 10;
% H = H(q)
E_q = zeros(2*N, B_max);
for b = 1:B_max
    H_q = [EA*eye(N), sqrt(b)*a; sqrt(b)*ad, EB*eye(N)];
    E_q(:, b) = eig(H_q);
end
fig_q = figure();
ax_q = axes(fig_q);

for i = 1:5
    plot(ax_q, 1:B_max, E_q(N+i,:));
    hold(ax_q, 'on');
end

% H = H(k)
% K_const = (10, 10)
Kx = 10;
Ky = 20;
E_k = zeros(2*N, B_max);
for b = 1:B_max
    H_k = [EA*eye(N), sqrt(b)*a- Kx + 1i*Ky; sqrt(b)*ad - Kx - 1i*Ky, EB*eye(N)];
    E_k(:, b) = eig(H_k);
end
fig_k = figure();
ax_k = axes(fig_k);
for i = 1:5
    plot(ax_k, 1:B_max, E_k(N+i,:));
    hold(ax_k, 'on');
end