%parpool('processes', 8);
tic;
w = 110e-3 * 1.6e-19;
theta = 2.5;
vF = 1e6;
eta = 1;
N_ll = 10;
p_max = 4;
q_max = 1;
p_range = 2:1:p_max;
q_range = 1:1:q_max;
coprime_pq_pairs = [];
for p = p_range
    for q = q_range
        if iscoprime([p, q]) && 2*p >= q
            coprime_pq_pairs = [coprime_pq_pairs, [p; q]];
        end
    end
end  

fig1 = figure();              
ax1 = axes(fig1);
data = [];
parfor ind = 1:length(coprime_pq_pairs)
    p = coprime_pq_pairs(1, ind);
    q = coprime_pq_pairs(2, ind);
    ll = LL(theta=theta, N_ll=N_ll, vF=vF, eta=eta, w=w);
    ll.precompute_init(p, q);
    ll.construct_H(ll.k1_max/2, ll.k2_max/2);
    cur_E = eig(ll.H)' / ll.e;
    ratio = 1/2 * q/p;
    cur_B = ll.B
    data = [data, [repelem(cur_B, length(cur_E)); cur_E]];
end
scatter(data(1,:), data(2,:), 5, 'k', 'filled');
toc;
