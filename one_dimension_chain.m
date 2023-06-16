function E = one_dimension_chain(unit_cell_size, k)
%ONE_DIMENSION_CHAIN Summary of this function goes here
%   Detailed explanation goes here
% this function calculates a 1-d chain of identical atoms
% with nn tb method.
% hopping term and lattice parameter are 1
% unit_cell_size is the number of atoms in a unit cell we consider,
% though they are all identical
% k input as momentum range
% return energy 
    l = numel(k);
    E = zeros(unit_cell_size, l);
    for ind = 1:l
        ck = k(ind);
        H = zeros(unit_cell_size);
        for i = 1:unit_cell_size
            if i - 1 <= 0
                H(i, i-1+unit_cell_size) = H(i, i-1+unit_cell_size) + exp(1i*ck);
            else
                H(i, i-1) = H(i, i-1) + exp(1i*ck);
            end
            if i + 1 > unit_cell_size
                H(i, i+1-unit_cell_size) = H(i, i+1-unit_cell_size) + exp(-1i*ck);
            else
                H(i, i+1) = H(i, i+1) + exp(-1i*ck);
            end
        end
    E(:, ind) = eig(H);
    end
end

