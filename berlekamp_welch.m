function msg = berlekamp_welch(R, num_err, mod_p)
%% FUNCTION
% This function uses the Berlekamp-Welch Algorithm to recover
% the original, uncorrupted message, MSG, of a received and potentially
% corrupted message vector R = r1 r2 r3 ...,
% which is expected to have at most NUM_ERR errors.
%
% MSG is returned by this function in the form
%                         m_1  m_2 ... m_n
% for a message of length n.
%
% Operations and messages are in modulo MOD_P, where MOD_P
% is prime.
%
%% More on Algorithm 
%
% The algorithm used by this function uses the following 
% equation to recover the original message, if recovery is possible:
%
% Let n be the length of the original message and
% k be the max number of errors (NUM_ERR) expected in the 
% corrupted message.
%           Q(i) = P(i)E(i) = (r_i)E(i)
% where 
%       * P is a polynomial of degree n - 1 representing 
%           the original message,
%       * E is a degree k polynomial whose roots are the 
%           erroneous elements of R, and
%       * Q is the resultant degree n + k - 1 polynomial
%           Q(x) =  [a_(n+k-1)]*x^(n+k-1) + .... 
%                   + [a_2]*x^2 + [a_1]*x + a_0
%
%       * r_i is the ith element of the received message R
% 
% From the equation above, n+2k linear equations can be generated to 
% solve for the coefficients of Q(x) and E(x).
%
% Then P(x) = Q(x)/E(x) since Q(x) = P(x)E(x) by construction.
% After finding P(x), which represents the original message,
% the original message can be recovered. This function assumes
% that P(x) is a coefficient representation of the original message.
%
% 
%% CODE 

% Check that R is a vector
if(~isvector(R))
    error('R must be a (1-D) vector.');
end

% R is made into a column vector, if not already one
if (size(R, 2) > 1)
    R = R';
end

r_length = numel(R);        %Length of R
n = r_length - 2*num_err;   %Length of original msg
q_in = (1:r_length)';       %array of consecutive numbers
deg_q = n + num_err - 1;    %Degree of Q(x)

%System of equations for Q(1)...Q(n+2k)
Q = [zeros(r_length, deg_q - 1), q_in, ones(r_length, 1)];
for i = (deg_q-1):-1:1
    Q(:, i) = Q(:, i+1).*q_in;
end

%System of equations for E(1)...E(n+2k)
E = [zeros(r_length, num_err-1), q_in, ones(r_length, 1)];
for i = num_err-1:-1:1
    E(:, i) = E(:, i+1).*q_in;
end

% r_iE(i)
for i = 1:(num_err+1)
    E(:,i) = R.*E(:,i);
end

% Q is set-up such that all unknowns a_i from Q(x) and b_i from E(x)
% are "on the LHS" so that the unknowns can be solved for via
% row-reduction.
Q = [Q, -1.*E(:,2:end), E(:,1)];
disp('Q, Before Mod:')
disp(Q);

Q = mod(Q, mod_p);
disp('Q, After Mod:')
disp(Q);

Q = row_reduce_p(Q, mod_p);
disp('Q, Row Reduced:')
disp(Q);

% z = [a_n ... a0 b_k-1 ... b0]
z = Q(:, end);
disp('z = [a_n ... a0 b_k-1 ... b0]:');
disp(z);

Q = z(1:deg_q+1); 
disp('Q(x):');
disp(Q);
E = [1; z(deg_q+2:end)];
disp('E(x):')
disp(E);

% P(x) = Q(x)/E(x)
[msg r] = deconv(Q,E);

% If there exists a remainder, there were more errors than
% buffered for or some another error occurred.
if (any(mod(r, 11) ~= 0))
    error('Message could not be recovered; more than k errors. ');
end

% Format msg from column vector to row vector of form m_1 m_2 ... m_n
msg = (flipud(mod(msg, mod_p)))';