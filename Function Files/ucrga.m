function R = ucrga(A)
% Unit consistent RGA
%{
----------------------------- Description ---------------------------------
RGA analysis where the matrix generalized inverse is the unit consistent
inverse, and not the Moore-Penrose inverse.

Copied directly from Uhlmann's paper - 
Uhlmann - 2019 - On the Relative Gain Array (RGA) with singular and
rectangular matrices

----------------------------- Inputs --------------------------------------
A : Input matrix, ideally this would be the plant gains (transfer function)
    at a certain frequency

----------------------------- Outputs -------------------------------------
R : RGA Matrix = A . ucinv(A)
    Here '.' indicates the Schur product, i.e, element by element product
------------------------------ Versions -----------------------------------
v1 : Copied into MATLAB by Suraj R Pawar, 6-5-2020
    - Initialize
%}

tol = 1e-15;
[m, n] = size(A);
L = zeros(m, n); M = ones(m, n);
S = sign(A); AA = abs(A);
idx = find(AA > 0.0); L(idx) = log(AA(idx));
idx = setdiff(1 : numel(AA), idx);
L(idx) = 0; M(idx) = 0;
r = sum(M, 2); c = sum(M, 1);
u = zeros(m, 1); v = zeros(1, n);
dx = 2*tol;
while (dx > tol)
idx = c > 0;
p = sum(L(:, idx), 1) ./ c(idx);
L(:, idx) = L(:, idx) - repmat(p, m, 1) .* M(:, idx);
v(idx) = v(idx) - p; dx = mean(abs(p));
idx = r > 0;
p = sum(L(idx, :), 2) ./ r(idx);
L(idx, :) = L(idx, :) - repmat(p, 1, n) .* M(idx, :);
u(idx) = u(idx) - p; dx = dx + mean(abs(p));
end
dl = exp(u); dr = exp(v);
S = S.* exp(L);
R = A .* transpose(pinv(S) .* (dl * dr)');
end


