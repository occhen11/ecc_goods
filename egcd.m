function [d, a, b] = egcd(x, y)

% Returns D, the GCD of X and Y, and A and B such that AX + BY = D

if (y == 0)
    d = x;
    a = 1;
    b = 0;
elseif (y > x)
    [d, b, a] = egcd(y, x);
else
    [d, a, b] = egcd(y, mod(x, y));
    temp = b;
    b = a - floor(x/y)*b;
    a = temp;
end
