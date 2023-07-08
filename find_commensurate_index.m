function [ans_m, ans_n] = find_commensurate_index(theta)
%FIND_COMMENSURATE_INDEX Summary of this function goes here
%   Detailed explanation goes here
% theta in degree unit
    % tolenrance in degree
    tol = 0.005;
    % store all possible m, n
    store = [];
    for n = 1:200
        for m = n:200
            cur_theta = acos((n^2 + 4*n*m + m^2)/(2*(n^2 + n*m + m^2)));
            cur_theta = cur_theta * 180 / pi;
            if abs(cur_theta - theta) <= tol
                store = [store, [m; n]];
                break;
            end
        end
    end
    if numel(store) == 0
        ME = MException('find_commensurate_index:notFound', 'Can not find a pair of indices satisfying the angle.');
        throw(ME);
    end
    [~, min_ind] = min(abs(store(:, 1) - store(:, 2)));
    ans_m = store(1, min_ind);
    ans_n = store(2, min_ind);
end

