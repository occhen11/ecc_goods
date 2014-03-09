function A = row_reduce_p(A, mod_p)

% This function row reduces matrix A to (reduced) row echelon form modulo p.

[r, c] = size(A);

for i = 1:r
    if(any(A(i,:)))
        for j = i+1:r
            A(j,:) = A(i,i).*A(j, :) - A(j,i).*A(i, :);
            A(j,:) = mod(A(j, :), mod_p);
        end
        [~, a] = egcd(A(i,i), mod_p);
        A(i,:) = mod(a*A(i, :), mod_p);
    end
end

for i = r:-1:1
    if(any(A(i,:)))
        for j = i-1:-1:1
           A(j,:) = A(i,i).*A(j, :) - A(j,i).*A(i, :);
        end
        A(i,:) = mod(A(i, :), mod_p);
    end
end

A(:, c) = mod(A(:, c), mod_p);