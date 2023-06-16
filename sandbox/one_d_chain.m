def 
unit_cell_size = 6;
H = zeros(unit_cell_size);
for i = 1:unit_cell_size
    if i - 1 <= 0
        H(i, i - 1 + unit_cell_size) = H(i, i - 1 + unit_cell_size) + 1 + 1i;
    else
        H(i, i - 1) = H(i, i - 1) + 1 + 1i;
    end
    if i + 1 > unit_cell_size
        H(i, i + 1 - unit_cell_size) = H(i, i + 1 - unit_cell_size) + 1 - 1i;
    else
        H(i, i + 1) = H(i, i + 1) + 1 - 1i;
    end
end
H
eig(H)