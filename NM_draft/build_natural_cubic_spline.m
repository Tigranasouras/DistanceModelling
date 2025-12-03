function coeffs = build_natural_cubic_spline(x, y)
% build_natural_cubic_spline
% Constructs a *natural* cubic spline through (x_i, y_i)
% following the 3.4 class notes notation.
%
% INPUT:
%   x : vector of knots  x_0 < x_1 < ... < x_n
%   y : vector of data   y_i = f(x_i)
%
% OUTPUT:
%   coeffs : (n x 4) matrix with rows [a_i, b_i, c_i, d_i]
%            defining
%            S_i(t) = a_i + b_i (t-x_i) + c_i (t-x_i)^2 + d_i (t-x_i)^3
%            on [x_i, x_{i+1}]

    x = x(:);  y = y(:);
    n = length(x) - 1;          % number of intervals

    % delta_i = x_{i+1} - x_i,   Delta_i = y_{i+1} - y_i
    delta  = diff(x);           % size n
    Delta  = diff(y);           % size n

    % --- Build (n+1)x(n+1) banded matrix A and RHS b for c_0..c_n ---
    A = zeros(n+1, n+1);
    b = zeros(n+1, 1);

    % Natural BC: c_0 = 0, c_n = 0
    A(1,1)       = 1;   b(1)     = 0;
    A(n+1,n+1)   = 1;   b(n+1)   = 0;

    % Interior equations: i = 1..n-1 (corresponding to nodes x_i)
    for i = 2:n
        A(i, i-1) = delta(i-1);
        A(i, i)   = 2*(delta(i-1) + delta(i));
        A(i, i+1) = delta(i);

        b(i) = 3*( Delta(i)/delta(i) - Delta(i-1)/delta(i-1) );
    end

    % Solve for c_i (these match the c_i in the notes)
    c_nodes = A \ b;  % size n+1, entries c_0..c_n

    % --- Compute a_i, b_i, c_i, d_i for each subinterval ---
    coeffs = zeros(n,4);  % [a_i, b_i, c_i, d_i]

    for i = 1:n
        a_i = y(i);
        c_i = c_nodes(i);          % our quadratic coefficient
        d_i = (c_nodes(i+1) - c_nodes(i)) / (3*delta(i));
        b_i = Delta(i)/delta(i) - (c_nodes(i+1) + 2*c_nodes(i))*delta(i)/3;

        coeffs(i,:) = [a_i, b_i, c_i, d_i];
    end
end
