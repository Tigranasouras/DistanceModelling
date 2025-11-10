% Program 0.1 Nested Multiplication
% Evaluates degree d polynomial from nested form using Horner's method with
% base points
% y = c_0 + (x-b_0)(c_1 + (x-b_1)(c_2 + (x-b_2)(c_3 ...
% Input: array of d+1 coefficients c (leading coefficient first),
%        x-coordinate x at which to evaluate,
%        array of d base points b (highest degree first)
% Output: value y of polynomial at x
function y = nest_base(c, x, b)

% base points are optional (assumed 0 if not provided; equiv. to nest.m if
% no base points are provided)
if nargin<3, b = zeros(length(c) - 1, 1); end

y = c(1) * ones(size(x));

for i = 2 : length(c)
  y = y .* (x - b(i - 1)) + c(i);
end

end
